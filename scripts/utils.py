#!/usr/bin/python
# -*- coding: utf-8 -*-
import csv, sys, re
def read_file(filename, sep ):
    with open(filename, "r") as load_file:
       reader = csv.reader(load_file, delimiter=sep)
       next(reader) # pour passer le header
       load_file=list(reader)
    return load_file
def read_file_withoutheader (filename, sep):
	with open(filename, "r") as load_file:
		reader = csv.reader(load_file, delimiter=sep)
		load_file=list(reader)
	return load_file
def create_output(table, columns):
	with open (table, 'w') as out:
		out.write(columns)
	return table