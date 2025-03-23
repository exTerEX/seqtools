from Bio import SeqRecord, SeqFeature


def extract_exact_feature(
    record: SeqRecord.SeqRecord,
    hit_start: int,
    hit_end: int,
    included_features: str | tuple = "CDS",
    verbose: bool = False,
) -> SeqFeature.SeqFeature:
    # For consistancy, this is converted to a tuple of strings
    if type(included_features) is str:
        included_features = (included_features,)

    if verbose:
        print(f"Include only these features: {included_features}")

    for feature in record.features:
        # Skip features not included in analysis
        if feature.type not in included_features:
            if verbose:
                try:
                    print(f"Skipped feature {feature.qualifiers['locus_tag']}")
                except KeyError:
                    print("Skipped feature without locus_id")
                    print(feature)

            continue

        feature_start = feature.location.start
        feature_end = feature.location.end

        # BUG: Non-standard genbank files may incorrectly use "gene" instead of "source"
        if hit_start >= feature_start and hit_end <= feature_end:
            if verbose:
                print(
                    f"Found feature in accession {feature.id}:{feature_start}..{feature_end}"
                )

            return feature
        else:
            if verbose:
                print("No feature in the given location")

            return SeqFeature.SeqFeature()  # UNTESTED: Check if this work as intended
