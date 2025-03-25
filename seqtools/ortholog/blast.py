import multiprocessing
from functools import partial

import numpy
import pandas

from seqtools.entrez import fetch_genbank
from seqtools.feature import extract_exact_feature
from seqtools.parser import ncbi_qaccver_header_split

# Constants

# The default data structure for tabular output from NCBI blastn.
DEFAULT_12COL_STRUCT = {
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
    "bitscore": numpy.dtype("float64"),
}

# The extended default data structure for tabular output from NCBI blastn.
DEFAULT_25COL_STRUCT = {  # FIXME: Update to correct data types
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
    "bitscore": numpy.dtype("float64"),
    "sallseqid": numpy.dtype("O"),
    "score": numpy.dtype("int64"),
    "nident": numpy.dtype("int64"),
    "positive": numpy.dtype("int64"),
    "gaps": numpy.dtype("int64"),
    "ppos": numpy.dtype("float64"),
    "qframe": numpy.dtype("O"),
    "sframe": numpy.dtype("O"),
    "qseq": numpy.dtype("O"),
    "sseq": numpy.dtype("O"),
    "qlen": numpy.dtype("int64"),
    "slen": numpy.dtype("int64"),
    "salltitles": numpy.dtype("O"),
}


def _iter_row_record(
    df: pandas.DataFrame,
    qaccver_col: str,
    saccver_col: str,
    sstart_col: str,
    send_col: str,
    update_query: bool,
) -> pandas.DataFrame:
    for index, row in df.iterrows():
        if update_query:
            qaccver, qstart, qend = ncbi_qaccver_header_split(row.loc[qaccver_col])

            # Fetch the genbank file assosiated with qaccver, and extract a feature
            # from the genbank file spanning the target region from qstart to qend.
            qfeature = extract_exact_feature(
                fetch_genbank(qaccver, cache=True), qstart, qend
            )

            # Update the dataframe with data data from qfeature
            df.at[index, "qtag"] = qfeature.qualifiers.get("locus_tag")
            df.at[index, "qstart"] = qfeature.location.start
            df.at[index, "qend"] = qfeature.location.end

        # Extract BLAST sequence hit information from the row in the dataframe.
        saccver = row.loc[saccver_col]
        sstart = row.loc[sstart_col]
        send = row.loc[send_col]

        # Fetch the genbank file assosiated with saccver.
        srecord = fetch_genbank(saccver, cache=True)

        # Extract a feature spanning the hit region from sstart to send.
        sfeature = extract_exact_feature(srecord, sstart, send)

        # Collect relevant data to the dataframe when sfeature is not NoneType,
        # Otherwise remove rows from the input dataframe.
        if sfeature.location is not None:
            df.at[index, "sstart"] = sfeature.location.start
            df.at[index, "send"] = sfeature.location.end
            df.at[index, "sframe"] = sfeature.location.strand
            df.at[index, "sseq"] = sfeature.location.extract(srecord)

            df.at[index, "stag"] = sfeature.qualifiers.get("locus_tag")[0]
            df.at[index, "psudo"] = sfeature.qualifiers.get("psudo")
        else:
            df.drop(index, inplace=True)  # UNTESTED: Unclear logic

    return df


def extract_gene_info(
    df: pandas.DataFrame,
    column_header: list | tuple | pandas.core.indexes.base.Index = None,
    n_threads: int = None,
    qaccver_col: str = None,  # TODO: Make handling alternative headers less complex
    saccver_col: str = None,
    sstart_col: str = None,
    send_col: str = None,
) -> pandas.DataFrame:
    # Use all available threads if n_threads is not specified.
    if not n_threads:
        n_threads = multiprocessing.cpu_count()
    # Use all available threads if user request more then available.
    elif multiprocessing.cpu_count() <= n_threads:
        n_threads = multiprocessing.cpu_count()

    # Check if standard blast output structure and standardize
    if list(df.dtypes) == list(DEFAULT_12COL_STRUCT.values()):
        df.columns = DEFAULT_12COL_STRUCT.keys()

        # Drop columns of little interest
        df = df.drop(["pident", "mismatch", "gapopen", "evalue", "bitscore"], axis=1)

        # Add columns necessary for analysis
        df["sframe"] = pandas.Series()
        df["sseq"] = pandas.Series()

        # Column headers
        qaccver_col = "qaccver"
        saccver_col = "saccver"
        sstart_col = "sstart"
        send_col = "send"
    elif list(df.dtypes) == list(DEFAULT_25COL_STRUCT.values()):
        df.columns = DEFAULT_25COL_STRUCT.keys()

        # Drop columns of little interest
        df = df.drop(
            [
                "pident",
                "mismatch",
                "gapopen",
                "evalue",
                "bitscore",
                "sallseqid",
                "score",
                "nident",
                "positive",
                "gaps",
                "ppos",
                "qframe",
                "qseq",
                "qlen",
                "slen",
                "salltitles",
            ],
            axis=1,
        )

        # Column headers
        qaccver_col = "qaccver"
        saccver_col = "saccver"
        sstart_col = "sstart"
        send_col = "send"
    else:
        df.columns = column_header
        # TODO: More checks

    # TODO: Check if custom col values are correct dtypes

    # Add columns for filling
    df["stag"] = pandas.Series()
    df["qtag"] = pandas.Series()
    df["psudo"] = pandas.Series()

    # Split dataframe in equally sized chunks
    chunks = numpy.array_split(df, n_threads)

    # Create a pool with 'n_threads' processes
    pool = multiprocessing.Pool(processes=n_threads)

    # Analyze each row in the dataframe by chunk in parallel using
    # the _iter_row_record function, with custom column names passed
    # as arguments secondary arguments.
    result = pool.map(
        partial(
            _iter_row_record,
            qaccver_col=qaccver_col,
            saccver_col=saccver_col,
            sstart_col=sstart_col,
            send_col=send_col,
            update_query=True,  # TODO: Give this as input
        ),
        chunks,
    )

    # Concatenate each result chunck into an updated dataframe
    df_result = pandas.concat(result)

    # Close pool for new requests
    pool.close()

    # Exit pool when all jobs are finished
    pool.join()

    return df_result
