#! /usr/bin/env python3

import argparse
import csv

parser = argparse.ArgumentParser(
    description="Find levenshtein distances between lists of strings."
)
parser.add_argument("strings", type=str, help="File containing list of strings")

args = parser.parse_args()


def edit_distance(x: str, y: str) -> int:
    """Wagner-Fischer Algorithm"""
    dist = [[0 for j in range(len(y) + 1)] for i in range(len(x) + 1)]

    for i in range(1, len(x) + 1):
        dist[i][0] = i

    for j in range(1, len(y) + 1):
        dist[0][j] = j

    for j in range(0, len(y)):
        for i in range(0, len(x)):
            cost = 1
            if x[i] == y[j]:
                cost = 0
            dist[i + 1][j + 1] = min(
                dist[i][j + 1] + 1, dist[i + 1][j] + 1, dist[i][j] + cost
            )
    return dist[len(x)][len(y)]


with open(args.strings, "r") as infile:
    string_list = [i[0].strip() for i in csv.reader(infile)]

min_dist = edit_distance(string_list[0], string_list[1])
for i in range(len(string_list)):
    for j in range(i + 1, len(string_list)):
        distance = edit_distance(string_list[i], string_list[j])
        if distance < min_dist:
            min_dist = distance

print(min_dist)
