import os
from itertools import combinations

############################################################################
# sORF extractor settings
############################################################################
MAX_SORF_LENGTH: int = 303  # 153
MAX_SPROT_LENGTH: int = (MAX_SORF_LENGTH - 3) // 3  # 50
MIN_SORF_LENGTH: int = 24  # https://legacy.uniprot.org/uniprot/?query=taxonomy%3A%22Bacteria+%5B2%5D%22+length%3A%5B*+TO+100%5D+fragment%3Ano+AND+reviewed%3Ayes&sort=length&desc=no
MIN_SPROT_LENGTH: int = (MIN_SORF_LENGTH - 3) // 3  # 7
UPSTREAM_LENGTH: int = 60
DOWNSTREAM_LENGTH: int = 150
RBS_UPSTREAM_LENGTH: int = 25  # 21 is the longest prodigal rbs_spacer + motif
START_CODONS: frozenset[str] = frozenset(('ATG', 'GTG', 'TTG', 'ATA', 'ATC', 'ATT', 'CTG'))
STOP_CODONS: frozenset[str] = frozenset(('TGA', 'TAG', 'TAA'))
TRANSLATION_TABLES: tuple[int, int, int] = (11, 4, 25)
############################################################################
# sORF product settings
############################################################################
HYPOTHETICAL_TERMS: frozenset[str] = frozenset(('hypothetical', 'putative', 'unknown', 'possible', 'uncharacterized', 'uncharacterised',
                                                'probable', 'dubious', 'doubtful', 'questionable'))
HYPOTHETICAL_PRODUCT_NAMES: set[str] = {'protein'}
for term in HYPOTHETICAL_TERMS:
    HYPOTHETICAL_PRODUCT_NAMES.add(term)
    HYPOTHETICAL_PRODUCT_NAMES.add(f'{term} protein')
for terms in combinations(HYPOTHETICAL_TERMS, 2):
    HYPOTHETICAL_PRODUCT_NAMES.add(f'{terms[0]} {terms[1]} protein')
HYPOTHETICAL_PRODUCT_NAMES: frozenset[str] = frozenset(HYPOTHETICAL_PRODUCT_NAMES)

############################################################################
# sORF homology settings
############################################################################
SRV_CUTOFF_VALID: float = 0.7  # hits with a srv below are discarded for validation
############################################################################
# sORF DB settings
############################################################################
SORFDB_VER: float = 1.0
HMM_MAX_SPROT_LENGTH: int = 50
############################################################################
# reference and representative genomes
############################################################################
REFERENCE_ACCESSIONS: frozenset[str] = frozenset(('GCA_000006945',   # Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
                                                  'GCA_000195955',   # Mycobacterium tuberculosis H37Rv
                                                  'GCA_000009045',   # Bacillus subtilis subsp. subtilis str. 168
                                                  'GCA_000005845',   # Escherichia coli str. K-12 substr. MG1655
                                                  'GCA_000008865',   # Escherichia coli O157:H7 str. Sakai
                                                  'GCA_000013425',   # Staphylococcus aureus subsp. aureus NCTC 8325
                                                  'GCA_000006765',   # Pseudomonas aeruginosa PAO1
                                                  'GCA_000240185',   # Klebsiella pneumoniae subsp. pneumoniae HS11286
                                                  'GCA_000196035',   # Listeria monocytogenes EGD-e
                                                  'GCA_000009085',   # Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819
                                                  'GCA_000006925',   # Shigella flexneri 2a str. 301
                                                  'GCA_000022005',   # Caulobacter vibrioides NA1000
                                                  'GCA_000191145',   # Acinetobacter pittii PHEA-2
                                                  'GCA_000008725',   # Chlamydia trachomatis D/UW-3/CX
                                                  'GCA_000007765'))  # Coxiella burnetii RSA 493

REPRESENTATIVE_ACCESSIONS: set[str] = set()
with open('{}/representative_genomes.txt'.format('/'.join(os.path.abspath(__file__).split('/')[:-1]))) as f:
    for line in f.readlines():
        REPRESENTATIVE_ACCESSIONS.add(line.strip())
REPRESENTATIVE_ACCESSIONS: frozenset[str] = frozenset(REPRESENTATIVE_ACCESSIONS)

# EOF
