import re
import logging
from pretextannotate.fetch_reports import fetch_sequence_reports

logger = logging.getLogger('pretextannotation_logger')

def custom_sort_order(molecule):
    """
    Sort chromosomes by numeric prefix then suffix, followed by named sex/autosome variants:
    - Pure digits (1,2,3...) in numeric order
    - Digits plus letter suffix (e.g. '1A','2B') by numeric then suffix
    - Named chromosomes: X, X1, X2, Y, W, Z, Z1, Z2, B, B1, B2
    - All others last
    """
    m = re.match(r"^(\d+)([A-Za-z]*)$", molecule)
    if m:
        return (int(m.group(1)), m.group(2))
    order_map = {
        'X': (1000, ''), 'X1': (1000, '1'), 'X2': (1000, '2'),
        'Y': (2000, ''),
        'W': (3000, ''),
        'Z': (4000, ''), 'Z1': (4000, '1'), 'Z2': (4000, '2'),
        'B': (5000, ''), 'B1': (5000, '1'), 'B2': (5000, '2')
    }
    return order_map.get(molecule, (float('inf'), molecule))


def extract_chromosomes_only(accession) -> list[dict]:
    """
    Fetch sequence reports, include only:
      - role='assembled-molecule' and assigned_molecule_location_type='Chromosome'
      - role='unlocalized-scaffold' (assigned to a chromosome)
    Exclude any role='unplaced-scaffold'.

    Sum lengths per chromosome name, record the primary accession and GC% for true chromosomes.
    Returns list of dicts: {'INSDC', 'molecule', 'length'(Mb), 'GC'}.
    """
    reports: list[dict] = fetch_sequence_reports(accession)
    chr_data = {}
    for rec in reports:
        name = rec.get('chr_name')
        role = rec.get('role')
        loc = rec.get('assigned_molecule_location_type', '')
        if not name:
            continue

        # include only true chromosomes and their unlocalized scaffolds
        if (role == 'assembled-molecule' and loc == 'Chromosome') or role == 'unlocalized-scaffold':
            length = rec.get('length', 0)
            entry = chr_data.setdefault(name, {'length': 0, 'INSDC': None, 'GC': None})
            entry['length'] += length

            # record accession/GC only for genuine chromosomes
            if role == 'assembled-molecule' and loc == 'Chromosome':
                entry['INSDC'] = rec.get('genbank_accession')
                entry['GC'] = rec.get('gc_percent')

    # build final list skipping names without a primary accession
    chrom_list = []
    for name, info in chr_data.items():
        # exclude organelle molecules
        if name.upper() in {'MT', 'PLTD'}:
            continue
        if not info['INSDC']:
            continue
        chrom_list.append({
            'INSDC': info['INSDC'],
            'molecule': name,
            'length': round(info['length'] / 1e6, 2),
            'GC': info['GC']
        })
    chrom_list.sort(key=lambda x: custom_sort_order(x['molecule']))
    return chrom_list


def get_chromosome_lengths(accession):
    """
    Total base-pair length of all true chromosomes + their unlocalized scaffolds.
    """
    return sum(int(c['length'] * 1e6) for c in extract_chromosomes_only(accession))
