with open('incf_arg_summary.tab', 'r') as f_in, open('output.txt', 'w') as f_out:
    # Read the header row
    header = f_in.readline().strip().split('\t')
    # Write the header row to the output file
    f_out.write('\t'.join(header) + '\n')
    # Iterate over the data rows
    for line in f_in:
        # Split the row into columns
        row = line.strip().split('\t')
        # Replace non-"." values with their corresponding column header
        for i in range(1, len(row)):
            if row[i] != '.':
                row[i] = header[i]
            else:
                row[i] = 'NA'
        # Write the modified row to the output file
        f_out.write('\t'.join(row) + '\n')