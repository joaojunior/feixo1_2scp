def read_file_format_or_library(name_file):
    f = open(name_file, 'r')
    quantity_costs = 0
    quantity_rows = 0
    quantity_cols = 0
    number_row = 0
    costs = []
    incidence_matrix = []
    line = f.readline()
    line = line.split()
    m = int(line[0])
    n = int(line[1])
    while quantity_costs < n:
        line = f.readline()
        line = line.split()
        for cost in line:
            costs.append(int(cost))
            quantity_costs += 1
    while number_row < m:
        col = 0
        row = [0] * n
        line = f.readline()
        line = line.split()
        quantity_cols = int(line[0])
        while col < quantity_cols:
            line = f.readline()
            line = line.split()
            for item in line:
                row[int(item)-1] = 1
                col += 1
        incidence_matrix.append(row)
        number_row += 1
    f.close()
    return costs, incidence_matrix