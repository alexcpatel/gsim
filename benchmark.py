#!/usr/bin/python

from functools import reduce
from timeit import default_timer as timer
from math import sqrt
import argparse
import os, os.path, sys
import subprocess

PROGRAMS = ['./gsim-bad/gsim-bad',
            './gsim-barneshut/gsim-barneshut',
            './gsim-lockfree/gsim-lockfree']

def run(args, program, threads, benchmark, seed, r):
  tracefile = ""
  command = ""
  if program == './gsim-bad/gsim-bad':
    tracefile = "%s/%s_%d_%d_%d_%d_%d_r%d.txt" % \
      (args.f, program.split('/')[2], threads, benchmark[1],
       args.b, args.s, seed, r)
    command = "%s %d %d %d %d %d &> %s" % \
      (program, threads, benchmark[1], args.b, args.s,
       seed, tracefile)
  else:
    tracefile = "%s/%s_%d_%d_%d_%d_%.2f_%d_r%d.txt" % \
      (args.f, program.split('/')[2], threads, benchmark[1],
       args.b, args.s, benchmark[0], seed, r)
    command = "%s %d %d %d %d %.2f %d &> %s" % \
      (program, threads, benchmark[1], args.b, args.s,
       benchmark[0], seed, tracefile)

  start = timer()
  if subprocess.call(command, shell=True):
    error = 'Encountered error while running %s.' % command
    raise Exception(error)
  return timer() - start

def benchmark(args, outfile):
  seed = int(timer())
  benchmarks = zip([val for val in [0.1,0.3,0.5] for _ in range(3)],
                   [args.b, int(sqrt(args.b)), 1]*3)

  start = timer()
  for program in PROGRAMS:
    for threads in range(1, args.t+1):
      for benchmark in benchmarks:
        result = reduce(min, map(lambda r: \
          run(args, program, threads, benchmark, seed, r), range(args.r)))
        resultstr = "BENCHMARK %14s: threads=%-2d clusters=%-5d bodies=%-5d "  \
          "steps=%-3d theta=%.2f seed=%d --- MS: %11.5f NPM: %11.5f\n" % \
          (program.split('/')[2], threads, benchmark[1], args.b, args.s,
           benchmark[0], seed, 1000.0 * result, (1.0e9 * result) / (args.b * args.s))
        outfile.write(resultstr)
        sys.stdout.write(resultstr);
      outfile.write('\n')
      sys.stdout.write('\n');
    outfile.write('\n')
    sys.stdout.write('\n');
  end = timer()
  endstr = "BENCHMARK TOTAL TIME: %f ms\n" % ((end-start) * 1000.0)
  outfile.write(endstr)
  sys.stdout.write(endstr);

def main():
  # parse command line arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', metavar='max threads', type=int, help='maximum thread count', default=16)
  parser.add_argument('-b', metavar='bodies', type=int, help='number of galaxy bodies', default=1000)
  parser.add_argument('-s', metavar='steps', type=int, help='number of simulation steps', default=100)
  parser.add_argument('-r', metavar='runs', type=int, help='number of runs', default=4)
  parser.add_argument('-f', metavar='tracepath', type=str, help='trace folder path', default='./traces')
  parser.add_argument('-o', metavar='outfile', type=str, help='output file path', default='./benchmark.txt')
  args = parser.parse_args()

  # initialize trace directory
  try:
    os.mkdir(args.f)
  except Exception:
    try:
      map(os.unlink, (os.path.join(args.f,f) for f in os.listdir(args.f)))
    except Exception as e:
      sys.stdout.write("%s\n" % e)
      sys.exit(-1)

  # open output file and run benchmarks
  try:
    with open(args.o, 'w') as outfile:
      benchmark(args, outfile)
  except Exception as e:
    sys.stdout.write("%s\n" % e)
    sys.exit(-1)

  sys.exit(0)

if __name__ == '__main__':
  main()
