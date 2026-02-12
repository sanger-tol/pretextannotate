import os
import json
import logging
from urllib import request as _request, error as _error, parse as _parse

logger = logging.getLogger('pretextannotation_logger')

# Setup Entrez credentials
entrez_email = os.getenv('ENTREZ_EMAIL', 'default_email')
entrez_api_key = os.getenv('ENTREZ_API_KEY', 'default_api_key')


def fetch_sequence_reports(accession):
    """
    Fetch all sequence reports from NCBI Datasets API and return assembled-molecule entries.
    urllib-only implementation with basic timeout and error handling.
    """
    base_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/sequence_reports"

    # Attach API key if available
    params = {}
    if entrez_api_key and entrez_api_key != "default_api_key":
        params["api_key"] = entrez_api_key
    url = base_url if not params else f"{base_url}?{_parse.urlencode(params)}"

    headers = {
        "accept": "application/json",
        "User-Agent": f"PretextAnnotate; {entrez_email}",
    }

    req = _request.Request(url, headers=headers, method="GET")

    try:
        with _request.urlopen(req, timeout=15) as resp:
            data_bytes = resp.read()
        data = json.loads(data_bytes.decode("utf-8", errors="replace"))
    except _error.HTTPError as e:
        status = getattr(e, "code", None)
        reason = getattr(e, "reason", "")
        logger.error(f"[NCBI - fetch_seq_reports] HTTP error {status} fetching sequence reports for {accession}: {reason}")
        return []
    except _error.URLError as e:
        logger.error(f"[NCBI - fetch_seq_reports] Network error fetching sequence reports for {accession}: {getattr(e, 'reason', e)}")
        return []
    except json.JSONDecodeError as e:
        logger.error(f"[NCBI - fetch_seq_reports] Failed to parse JSON for {accession}: {e}")
        return []

    reports = data.get("reports", [])
    if not isinstance(reports, list):
        logger.error("[NCBI - fetch_seq_reports] Unexpected JSON structure: 'reports' is not a list")
        return []

    logger.error("[NCBI - fetch_seq_reports] Successfully called API")

    # include both assembled chromosomes and their unlocalized scaffolds
    return [r for r in reports if r.get("role") in ("assembled-molecule", "unlocalized-scaffold")]
