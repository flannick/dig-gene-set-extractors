# Paired pull requests

A wrapper submission records the DIG repository URL, exact commit, and the
workflow/converter identifiers it calls. Before merging paired changes, run
`geneset-extractors submission validate <identifier>` at that exact DIG commit.
The DIG contract is backwards compatible: it adds a new `submission` CLI
namespace and does not rename existing converter or workflow commands.
