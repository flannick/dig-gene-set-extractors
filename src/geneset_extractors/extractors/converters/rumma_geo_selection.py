"""Reproducibly select RummaGEO query records before membership reconstruction."""
from __future__ import annotations
import csv, json, re
from pathlib import Path
from geneset_extractors.core.metadata import input_file_record, make_metadata, write_metadata
from geneset_extractors.core.provenance import activate_runtime_context

CONTROL_TERMS = frozenset(('wt wildtype control cntrl ctrl uninfected normal untreated unstimulated shctrl ctl healthy sictrl sicontrol ctr wild dmso'.split()))
GENE_TERMS = ('knockout', 'crispr ko', 'overexpression', 'inhibition', 'knockdown')
DRUG_EXCLUSIONS = frozenset(('1B ATPA AVA C-1 CDC FIT ITE PIT RITA TRIM compe iq niacin pen rutin'.split()))

def is_control(value: str) -> bool:
    return any(re.sub('[^a-zA-Z]', '', word) in CONTROL_TERMS for word in value.lower().split())

def _nodes(path: Path):
    data=json.loads(path.read_text(encoding='utf-8'))
    return data['records'] if isinstance(data,dict) and 'records' in data else data

def run(args):
    activate_runtime_context('rumma_geo_selection', getattr(args,'provenance_overlay_json',None))
    query, drugs = Path(args.query_records_json).resolve(), Path(args.drug_terms_json).resolve() if args.drug_terms_json else None
    records=_nodes(query); model=args.model_id
    if model == 'HZ1' and drugs is None: raise ValueError('HZ1 requires --drug_terms_json')
    terms=set(GENE_TERMS) if model=='HZ2' else set(json.loads(drugs.read_text(encoding='utf-8')))
    selected={}
    for r in sorted(records, key=lambda x:(str(x.get('search_term','')),str(x.get('id','')))):
        if r.get('search_term') not in terms or (model=='HZ1' and r.get('search_term') in DRUG_EXCLUSIONS): continue
        try:
            source,direction=str(r['term']).split(); gse,c1,_,c2,species=source.split('-'); titles=r['sampleGroups']['titles']; a,b=titles[c1],titles[c2]
        except (KeyError,ValueError): continue
        status='signature'
        if is_control(a) and is_control(b): status='2 control'
        elif not is_control(a) and not is_control(b): status='2 pert'
        elif is_control(a): status='reversed'
        if status not in ('signature','reversed'): continue
        uid=str(r['id'])
        selected.setdefault(uid, {'uuid':uid,'source_term':source,'model_id':model,'status':status,'search_term':r['search_term'],'gse':gse,'condition_1':c1,'condition_2':c2,'species':species,'source_direction':direction,'condition_1_title':a,'condition_2_title':b,'title':str(r.get('title',''))})
    out=Path(args.out_dir).resolve(); out.mkdir(parents=True,exist_ok=True); rows=sorted(selected.values(),key=lambda r:(r['source_term'],r['uuid']))
    fields=list(rows[0]) if rows else ['uuid','source_term','model_id','status','search_term','gse','condition_1','condition_2','species','source_direction','condition_1_title','condition_2_title','title']
    with (out/'selection_manifest.tsv').open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
    (out/'query_records.used.json').write_text(json.dumps(rows,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    files=[input_file_record(query,'cached_rummageo_graphql_records')]+([input_file_record(drugs,'sigcom_lincs_drug_terms')] if drugs else [])
    meta=make_metadata('rumma_geo_selection',{'model_id':model,'control_terms':sorted(CONTROL_TERMS),'gene_search_terms':list(GENE_TERMS),'drug_exclusions':sorted(DRUG_EXCLUSIONS),'uuid_deduplication':'stable_first_by_search_term_then_uuid'},'metadata','rna_seq','human','hg38',files,{'mode':'none'},{'weight_type':'unweighted','normalization':{},'aggregation':'selection'}, {'n_input_features':len(records),'n_genes':0,'n_features_assigned':len(rows),'fraction_features_assigned':len(rows)/len(records) if records else 0,'n_gene_sets':len(rows)},output_files=[{'path':'selection_manifest.tsv','role':'selection_manifest'},{'path':'query_records.used.json','role':'selected_query_records'}],gene_set_description='RummaGEO cached-query selection manifest')
    write_metadata(out/'geneset.meta.json',meta); return {'n_peaks':len(rows),'n_genes':0,'out_dir':str(out)}
