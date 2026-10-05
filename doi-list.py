import ccdc.io
import ccdc.search
from bsubpc_match import TEMPLATE_SMARTS

template_substructure = ccdc.search.SMARTSSubstructure(TEMPLATE_SMARTS)

# Create the BsubPc substructure search
subpc_search = ccdc.search.SubstructureSearch()
subpc_search.add_substructure(template_substructure)

# Run the CSD search
hits = subpc_search.search()
dois = set()
for hit in hits:
    entry = hit.entry
    publication = entry.publication
    doi = publication.doi
    if doi is None:
        continue
    if doi not in dois:
        print(doi)
        dois.add(doi)
