from pathlib import Path

from Bio import Entrez, SeqIO, SeqRecord


def fetch_genbank(
    accession: str,
    cache: bool = False,
    local_only: bool = False,
    accession_cache_base_path: str = ".cache/entrez/",
    verbose: bool = False,
) -> SeqRecord.SeqRecord:
    # Create a path for loading or caching genbank files
    accession_cache_path = Path(accession_cache_base_path + accession + ".gbff")

    if accession_cache_path.exists():
        record = SeqIO.read(accession_cache_path, "genbank")

        if verbose:
            print(
                f"Loading genbank file for accession {accession} from local cache located in {accession_cache_path}"
            )
    elif accession_cache_path.exists() and local_only:
        record = SeqIO.read(accession_cache_path, "genbank")
    else:
        with Entrez.efetch(
            db="nuccore", id=accession, rettype="gbwithparts", retmode="text"
        ) as handle:
            record = SeqIO.read(handle, "genbank")

            if cache:
                accession_cache_path.parent.mkdir(parents=True, exist_ok=True)

                SeqIO.write(record, accession_cache_path, "genbank")

            if verbose and cache:
                print(
                    f"Fetching genbank file of {accession} from NCBI and storing a local copy in {accession_cache_path}"
                )
            elif verbose:
                print(f"Fetching genbank file of {accession} from NCBI")

    return record
