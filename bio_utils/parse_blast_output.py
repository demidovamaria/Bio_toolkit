def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parses a BLAST output TXT file to extract the top-match protein name (Description)
    from the first line of the alignments table for each Query, collects unique names,
    sorts them alphabetically, and writes to the output file.

    Arguments:
        input_file: Path to the input BLAST TXT file.
        output_file: Path to the output file where sorted protein names will be written.
    """
    top_protein_names = set()

    file = open(input_file, 'r')
    full_file = file.read()
    queries = full_file.split("Query #")

    for query in queries:
        next_line_best = False
        name_start_index = -1
        for line in query.splitlines():
            line = line.strip()
            if next_line_best:
                print('hit: ' + line)
                if name_start_index != -1:
                    top_protein_names.add(line[0:name_start_index])
                break
            if not line:  # skip empty
                continue
            if line.startswith("Description"):
                next_line_best = True
                name_start_index = line.find("Name")

        # Sort protein names alphabetically
    sorted_protein_names = sorted(top_protein_names)

    # Write to output file
    output = open(output_file, 'w')
    for name in sorted_protein_names:
        output.write(name + '\n')
