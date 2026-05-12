from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


Json = dict[str, Any]


@dataclass(frozen=True)
class Checksum:
    algorithm: str
    encoding: str
    value: str

    def to_jsonld(self) -> Json:
        return {
            "algorithm": self.algorithm,
            "encoding": self.encoding,
            "value": self.value,
        }


@dataclass
class FileEntity:
    id: str
    name: str
    file_role: str
    reproducibility_role: str
    entity_types: list[str]
    relative_path: str | None = None
    absolute_path_observed: str | None = None
    content_size: int | None = None
    checksums: list[Checksum] = field(default_factory=list)
    was_generated_by: str | None = None
    was_derived_from: str | None = None
    c2m2: Json | None = None
    extra: Json = field(default_factory=dict)

    def to_jsonld(self) -> Json:
        payload: Json = {
            "@id": self.id,
            "@type": self.entity_types,
            "name": self.name,
            "dig:file_role": self.file_role,
            "dig:reproducibility_role": self.reproducibility_role,
        }
        if self.relative_path:
            payload["dig:relative_path"] = self.relative_path
        if self.absolute_path_observed:
            payload["dig:absolute_path_observed"] = self.absolute_path_observed
        if self.content_size is not None:
            payload["contentSize"] = self.content_size
        if self.checksums:
            payload["dig:checksums"] = [item.to_jsonld() for item in self.checksums]
        if self.was_generated_by:
            payload["prov:wasGeneratedBy"] = {"@id": self.was_generated_by}
        if self.was_derived_from:
            payload["prov:wasDerivedFrom"] = {"@id": self.was_derived_from}
        if self.c2m2 is not None:
            payload["dig:c2m2"] = self.c2m2
        payload.update(self.extra)
        return payload


@dataclass
class CommandStep:
    id: str
    name: str
    step_index: int
    entrypoint: str
    observed_command: str
    replay_command: str
    command_argv: list[str]
    module: str | None = None
    working_directory: str = "."
    used: list[str] = field(default_factory=list)
    generated: list[str] = field(default_factory=list)
    parameter_set_id: str | None = None
    associated_with: list[str] = field(default_factory=list)
    extra: Json = field(default_factory=dict)

    def to_jsonld(self) -> Json:
        payload: Json = {
            "@id": self.id,
            "@type": ["prov:Activity", "dig:CommandStep", "dig:AnalysisStep"],
            "name": self.name,
            "dig:step_index": self.step_index,
            "dig:entrypoint": self.entrypoint,
            "dig:observed_command": self.observed_command,
            "dig:replay_command": self.replay_command,
            "dig:command_argv": self.command_argv,
            "dig:working_directory": self.working_directory,
            "prov:used": [{"@id": item} for item in self.used],
            "prov:generated": [{"@id": item} for item in self.generated],
            "prov:wasAssociatedWith": [{"@id": item} for item in self.associated_with],
        }
        if self.module:
            payload["dig:module"] = self.module
        if self.parameter_set_id:
            payload["dig:parameters"] = {"@id": self.parameter_set_id}
        payload.update(self.extra)
        return payload


@dataclass
class ResolvedParameterSet:
    id: str
    step_id: str
    cli_parameters: Json
    resolved_parameters: Json
    parameter_fingerprint: str
    parameter_overrides: list[Json] = field(default_factory=list)
    parameter_warnings: list[str] = field(default_factory=list)

    def to_jsonld(self) -> Json:
        return {
            "@id": self.id,
            "@type": ["prov:Entity", "dig:ResolvedParameterSet"],
            "dig:step": {"@id": self.step_id},
            "dig:cli_parameters": self.cli_parameters,
            "dig:resolved_parameters": self.resolved_parameters,
            "dig:parameter_fingerprint": self.parameter_fingerprint,
            "dig:parameter_overrides": self.parameter_overrides,
            "dig:parameter_warnings": self.parameter_warnings,
        }


@dataclass
class SoftwareRepository:
    id: str
    repo_url: str
    git_commit: str
    git_ref: str | None = None
    git_dirty: bool | None = None
    install_command: str | None = None

    def to_jsonld(self) -> Json:
        payload: Json = {
            "@id": self.id,
            "@type": ["prov:Entity", "dig:SoftwareRepository"],
            "dig:repo_url": self.repo_url,
            "dig:git_commit": self.git_commit,
        }
        if self.git_ref is not None:
            payload["dig:git_ref"] = self.git_ref
        if self.git_dirty is not None:
            payload["dig:git_dirty"] = self.git_dirty
        if self.install_command:
            payload["dig:install_command"] = self.install_command
        return payload


@dataclass
class SoftwareEnvironment:
    id: str
    platform: str
    python_version: str
    r_version: str | None = None
    container_image: str | None = None
    container_digest: str | None = None
    pip_freeze_file_id: str | None = None
    r_session_info_file_id: str | None = None
    r_packages: list[dict[str, str]] = field(default_factory=list)

    def to_jsonld(self) -> Json:
        payload: Json = {
            "@id": self.id,
            "@type": ["prov:Entity", "dig:SoftwareEnvironment"],
            "platform": self.platform,
            "dig:python_version": self.python_version,
            "dig:r_packages": self.r_packages,
        }
        if self.r_version:
            payload["dig:r_version"] = self.r_version
        if self.container_image:
            payload["dig:container_image"] = self.container_image
        if self.container_digest:
            payload["dig:container_digest"] = self.container_digest
        if self.pip_freeze_file_id:
            payload["dig:pip_freeze_file"] = {"@id": self.pip_freeze_file_id}
        if self.r_session_info_file_id:
            payload["dig:r_session_info_file"] = {"@id": self.r_session_info_file_id}
        return payload


@dataclass
class ReplayPlan:
    id: str
    scope: str
    entrypoint: str
    steps: list[str]
    required_inputs: list[str]
    expected_outputs: list[str]
    comparison_targets: list[str]

    def to_jsonld(self) -> Json:
        return {
            "@id": self.id,
            "@type": ["prov:Plan", "dig:ReplayPlan"],
            "dig:scope": self.scope,
            "dig:entrypoint": self.entrypoint,
            "dig:steps": [{"@id": item} for item in self.steps],
            "dig:required_inputs": [{"@id": item} for item in self.required_inputs],
            "dig:expected_outputs": [{"@id": item} for item in self.expected_outputs],
            "dig:comparison_targets": [{"@id": item} for item in self.comparison_targets],
        }


@dataclass
class GeneSetProvenanceRecord:
    record_id: str
    schema_version: str
    geneset: Json
    model: Json
    replay_plan_id: str
    graph: list[Any]

