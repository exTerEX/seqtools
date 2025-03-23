def find_upstream_bounderies(
    strand: int,
    start: int,
    stop: int,
    length: int = 200,
    include_first_codon: bool = True,
    verbose: bool = False,
) -> tuple[int, int]:
    if strand not in [-1, 1]:
        raise ValueError("Strand can only be -1 or 1")

    if strand == 1:
        start_region, stop_region = start - length - 1, start

        if verbose:
            print(f"Found region in sense from {start_region} to {stop_region}")

        if include_first_codon:
            stop_region += 2
    else:
        start_region, stop_region = stop, stop + length + 1

        if verbose:
            print(f"Found region in antisense from {start_region} to {stop_region}")

        if include_first_codon:
            start_region -= 3

    if include_first_codon and verbose:
        print("Included the first codon / 3 nucleotides of ")

    return int(start_region), int(stop_region)
