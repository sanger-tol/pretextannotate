#!/usr/bin/env python3


import os
import re
import requests
from Bio import Entrez

# Setup Entrez credentials
entrez_email = os.getenv('ENTREZ_EMAIL', 'default_email')
entrez_api_key = os.getenv('ENTREZ_API_KEY', 'default_api_key')


def fetch_sequence_reports(accession):
    """
    Fetch all sequence reports from NCBI Datasets API and return assembled-molecule entries.
    """
    api_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/sequence_reports"
    headers = {
        'accept': 'application/json',
        'User-Agent': f'PretextAnnotate; {entrez_email}'
    }
    resp = requests.get(api_url, headers=headers)
    resp.raise_for_status()
    reports = resp.json().get('reports', [])
    # include both assembled chromosomes and their unlocalized scaffolds
    return [r for r in reports if r.get('role') in ('assembled-molecule','unlocalized-scaffold')]


def get_longest_scaffold(accession):
    """
    Return the length (Mb) of the longest assembled-molecule (scaffold or chromosome).
    """
    reports = fetch_sequence_reports(accession)
    if not reports:
        return None
    longest = max(reports, key=lambda r: r.get('length', 0))
    return round(longest.get('length', 0) / 1e6, 2)


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
    reports = fetch_sequence_reports(accession)
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


def prim_chromosome_table(accession):
    """
    Return a sorted list of true chromosomes for primary assembly.
    """
    return extract_chromosomes_only(accession)


def combine_haplotype_chromosome_tables(hap1_acc, hap2_acc):
    """
    Return side-by-side rows for two haplotype chromosome lists.
    """
    hap1 = extract_chromosomes_only(hap1_acc)
    hap2 = extract_chromosomes_only(hap2_acc)
    hap1.sort(key=lambda x: custom_sort_order(x['molecule']))
    hap2.sort(key=lambda x: custom_sort_order(x['molecule']))
    combined = []
    for i in range(max(len(hap1), len(hap2))):
        combined.append({
            'hap1_INSDC': hap1[i]['INSDC'] if i < len(hap1) else "",
            'hap1_molecule': hap1[i]['molecule'] if i < len(hap1) else "",
            'hap1_length': hap1[i]['length'] if i < len(hap1) else "",
            'hap1_GC': hap1[i]['GC'] if i < len(hap1) else "",
            'hap2_INSDC': hap2[i]['INSDC'] if i < len(hap2) else "",
            'hap2_molecule': hap2[i]['molecule'] if i < len(hap2) else "",
            'hap2_length': hap2[i]['length'] if i < len(hap2) else "",
            'hap2_GC': hap2[i]['GC'] if i < len(hap2) else ""
        })
    return combined


def get_chromosome_lengths(accession):
    """
    Total base-pair length of all true chromosomes + their unlocalized scaffolds.
    """
    chroms = extract_chromosomes_only(accession)
    # extract_chromosomes_only lengths are in Mb
    return sum(int(c['length'] * 1e6) for c in chroms)


def calculate_percentage_assembled(info):
    """
    Compute percentage of total genome sequence assigned to chromosomes.
    Works for prim_alt and hap_asm.
    """
    asm = info.get('assemblies_type')
    if asm == 'prim_alt':
        total_chr = get_chromosome_lengths(info['prim_accession'])
        genome = float(info.get('genome_length_unrounded', 0))
        pct = round((total_chr / genome) * 100, 2) if genome else 0
        return {'total_chromosome_length': total_chr,
                'genome_length_unrounded': genome,
                'perc_assembled': pct}
    elif asm == 'hap_asm':
        total1 = get_chromosome_lengths(info['hap1_accession'])
        total2 = get_chromosome_lengths(info['hap2_accession'])
        g1 = float(info.get('hap1_genome_length_unrounded', 0))
        g2 = float(info.get('hap2_genome_length_unrounded', 0))
        p1 = round((total1 / g1) * 100, 2) if g1 else 0
        p2 = round((total2 / g2) * 100, 2) if g2 else 0
        return {'hap1_chromosome_length': total1,
                'hap2_chromosome_length': total2,
                'hap1_perc_assembled': p1,
                'hap2_perc_assembled': p2}
    else:
        raise ValueError("Unsupported assembly type")


def identify_sex_chromosomes(chr_list):
    """
    Return sorted list of sex chromosomes and their variants found.
    """
    valid = {'X','X1','X2','Y','W','Z','Z1','Z2','B','B1','B2', 'U', 'V'}
    found = {c['molecule'].upper() for c in chr_list if c['molecule'].upper() in valid}
    return sorted(found)
