#!/usr/bin/env python3

import argparse


def getargs():
    parser = argparse.ArgumentParser(
        description="Script safely downsample logscaled data."
    )

    options = parser.add_argument_group("Options")
    options.add_argument(
        "-i", "--input-csv", required=True, help="CSV file to be downsampled"
    )
    options.add_argument(
        "-n",
        "--max-points",
        required=True,
        type=str,
        help="Maximum number of points ot be included on the downsampled data",
    )
    options.add_argument(
        "-x",
        "--x-column",
        required=False,
        type=int,
        default=0,
        help="The index of the column with the x values to use for the distance calculations",
    )
    options.add_argument(
        "-y",
        "--y-column",
        required=False,
        default=1,
        type=int,
        help="The index of the column with the y values to use for the distance calculations",
    )

    return parser.parse_args()


def readinpoints(inputcvs):
    xs = []
    ys = []

    F = open(inputcvs, "r")
    header = F.readline()

    data = F.readline()
    while data:
        data = data.strip()
        cols = data.split(",")
        xs.append(float(cols[args.x_column]))
        ys.append(float(cols[args.y_column]))
        data = F.readline()
    F.close()
    return (xs, ys, header)


def finddist(x1, y1, x2, y2):
    dist = ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5

    return dist


def findmaxdist(xs, ys, N):
    minx = min(xs)
    maxx = max(xs)
    miny = min(ys)
    maxy = max(ys)
    totaldist = finddist(minx, miny, maxx, maxy)
    maxdist = totaldist / float(N)

    return maxdist


def downsample(xs, ys, maxdist, inputcsv, header):
    # Always keep the first point/row
    downxs = [xs[0]]
    downys = [ys[0]]
    accepted_indices = [0]  # A list of the indices of rows to keep

    for i in range(1, len(xs) - 1):
        mydist = finddist(downxs[-1], downys[-1], xs[i], ys[i])
        if mydist >= maxdist:
            downxs.append(xs[i])
            downys.append(ys[i])
            # Add this index to the list, we will keep this row
            accepted_indices.append(i)

    # Always add the last point/row
    downxs.append(xs[-1])
    downys.append(ys[-1])
    accepted_indices += [len(xs) - 1]

    # rounded version of downys for error checking
    rounded_downys = [round(yval, 4) for yval in downys]

    outfile = inputcsv.replace(".csv", "_downsampled.csv")
    with open(inputcsv, "r") as input_rows:
        with open(outfile, "w") as output_file:
            output_file.write(header)
            # skip the header in the input rows
            next(input_rows)
            i = 0
            # Note: i starts at 0 becuase the indices in accepted_indices correspond to the indices of the
            # xs and ys objects, which skip the header and use index 0 for the first row of data
            rows_added = 0
            for line in input_rows:
                if i in accepted_indices:
                    data = line.strip()
                    cols = data.split(",")
                    output_file.write(line)
                    # This probably adds a lot of time but I'm worried about off-by-one error
                    # Comparing an ys and downys values are floats, rounding in case there's a minor difference
                    # in one of the last decimal places
                    check_val = round(float(cols[args.y_column]), 4)
                    if not (check_val in rounded_downys):
                        raise AssertionError(
                            "off by one check failed.  col val is: "
                            + str(check_val)
                            + "\nexpected rounded downys value is: "
                            + str(rounded_downys[rows_added])
                            + "\ncols is: "
                            + ",".join(cols)
                            + "\ni is "
                            + str(i)
                        )

                i += 1
                rows_added += 1

    # F=open(outfile,"w")
    # F.write(header)

    # for i in range(len(downxs)):
    # 	F.write("{},{}\n".format(downxs[i],downys[i]))

    # F.close()


args = getargs()
xpoints, ypoints, header = readinpoints(args.input_csv)
maxdist = findmaxdist(xpoints, ypoints, args.max_points)
downsample(xpoints, ypoints, maxdist, args.input_csv, header)
