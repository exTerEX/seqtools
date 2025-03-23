# TODO: Check if input header contain all information necessary
def ncbi_qaccver_header_split(header: str) -> tuple[str, int, int]:
    # TODO: add try-expect / regex to catch differences between "-" and ".."
    accession, location = header.split(":")
    start, end = location.split("-")

    return str(accession), int(start), int(end)
