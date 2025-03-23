import multiprocessing

import numpy
import pandas

from seqtools.entrez import fetch_genbank
from seqtools.feature import extract_exact_feature
from seqtools.parser import ncbi_qaccver_header_split

# Constants
DEFAULT_12COL_STRUCT = {  # FIXME: Update to correct data types
    "qaccver": numpy.dtype("O"),
    "saccver": numpy.dtype("O"),
    "pident": numpy.dtype("float64"),
    "length": numpy.dtype("int64"),
    "mismatch": numpy.dtype("int64"),
    "gapopen": numpy.dtype("int64"),
    "qstart": numpy.dtype("int64"),
    "qend": numpy.dtype("int64"),
    "sstart": numpy.dtype("int64"),
    "send": numpy.dtype("int64"),
    "evalue": numpy.dtype("float64"),
    "bitscore": numpy.dtype("int64"),
}

# TODO: DEFAULT_25COL_STRUCT


def _iter_row_record(
    df: pandas.DataFrame,
) -> pandas.DataFrame:
    for index, row in df.iterrows():
        qaccver, qstart, qend = ncbi_qaccver_header_split(row.loc["qaccver"])

        # TODO: cache if header didn't contain a proper accession and location
        qrecord = fetch_genbank(qaccver, cache=True)

        qfeature = extract_exact_feature(qrecord, qstart, qend)

        # TODO: Write result back to dataframe using .set_value()

        # Extract sequence hit information
        saccver = row.loc["saccver"]
        sstart = row.loc["sstart"]
        send = row.loc["send"]

        # Fetch the correct genbank
        srecord = fetch_genbank(saccver, cache=True)

        sfeature = extract_exact_feature(srecord, sstart, send)

        # TODO: Write result back to dataframe using .set_value()

    return df


def extract_gene_info(
    df: pandas.DataFrame,
    columns: list | tuple | pandas.core.indexes.base.Index = None,
    qaccver_col: str = None,
    saccver_col: str = None,
    qstart_col: str = None,
    qend_col: str = None,
    sstart_col: str = None,
    send_col: str = None,
    n_threads: int = None,
) -> pandas.DataFrame:
    if not n_threads:  # Use all if not specified
        n_threads = multiprocessing.cpu_count()
    elif multiprocessing.cpu_count() <= n_threads:  # Less cores then requested
        n_threads = multiprocessing.cpu_count()

    # Check if standard blast output structure and standardize
    if list(df.dtypes) == list(DEFAULT_12COL_STRUCT.values()):
        df.columns = DEFAULT_12COL_STRUCT.keys()
    # TODO: 25 column dataframe
    # TODO: If non standard, check that "qaccver", "saccver", "sstart", "ssend" exist

    # Split dataframe in equally sized chunks
    chunks = numpy.array_split(df, n_threads)

    # Create a pool with 'n_threads' processes
    pool = multiprocessing.Pool(processes=n_threads)

    # Analyze each row in the dataframe by chunk in parallel
    result = pool.map(_iter_row_record, chunks)

    # Concatenate each result chunck into an updated dataframe
    df = pandas.concat(result)

    # Close pool for new requests
    pool.close()

    # Exit pool when all jobs are finished
    pool.join()

    return df
