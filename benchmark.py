#!/usr/bin/python

import argparse
import random
import datetime
import math
import os
import os.path
import sys
import subprocess
import time

seed = int(time.time())
random.seed(seed)

def outmsg(s, outFile, noreturn = False):
    if type(s) != str:
        s = str(s)
    if len(s) > 0 and s[-1] != '\n' and not noreturn:
        s += "\n"
    sys.stdout.write(s)
    sys.stdout.flush()
    if outFile is not None:
        outFile.write(s)

def doRun(cmd_args, benchmark, outFile):
    cmdline = " ".join(cmd_args)
    tstart = datetime.datetime.now()

    try:
        print("Running '%s'" % cmdline)
        simProcess = subprocess.Popen(cmd_args)
        simProcess.wait()
        returnCode = simProcess.returncode
#        stream = simProcess.communicate()
    except Exception as e:
        print("Execution of command '%s' failed with exception: %s" % (cmdline, e))
        return None
    if returnCode == 0:
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        return secs
    else:
        print("Execution of command '%s' failed with code: %s" % (cmdline, returnCode))
        return None

def bestRun(cmd_args, benchmark, runs, outFile):
    sofar = sys.maxsize
    for r in range(runs):
        if runs > 1:
            print("Run #%d" % (r+1))
        secs = doRun(cmd_args, benchmark, outFile)
        if secs == None:
            return None
        sofar = min(sofar, secs)
    return sofar

def run_benchmark(quick, clusters, bodies, N, runs, benchmark, seq_dir, mpi_dir, processes, outFile, theta):
    if processes > 1:
        cmd_args = [] # TODO
    else:
        cmd_args = [seq_dir, str(clusters), str(bodies), str(N), str(theta), str(seed)]
    secs = bestRun(cmd_args, benchmark, runs, outFile)

    if not secs == None:
        bmoves = bodies * N # TODO
        npm = 1e9 * secs/bmoves
        print(secs, npm)
        return secs, npm
    else:
        return None

def get_benchmark(benchmark, bodies):
    if benchmark == '1': # uniform
        clusters = bodies
        return clusters, bodies
    elif benchmark == '2': # one clump
        clusters = 1
        return clusters, bodies
    elif benchmark == '3': # moderate clumping
        clusters = bodies // 4
        return clusters, bodies

def run_benchmarks(quick, bodies, N, runs, benchmarks, seq_dir, mpi_dir, processes, outFile):
    for theta in [0, 0.5, 1.0]:
        outmsg("PARAMETERS ----- bodies: %d, steps: %d, runs: %d, processes: %d, theta: %f" %
                (bodies, N, runs, processes, theta), outFile)
        for benchmark in list(benchmarks):
            clusters, bodies = get_benchmark(benchmark, bodies)
            results = run_benchmark(quick, clusters, bodies, N, runs, benchmark, seq_dir, mpi_dir, processes, outFile, theta)
            outmsg("Benchmark %s -- Seconds: %f, NPM: %f" % (benchmark, results[0], results[1]), outFile)

def run(args):
    fname = args.f
    try:
        outFile = open(fname, 'w')
    except Exception as e:
        outFile = None
        outmsg("Couldn't open output file '%s' with exception %s" % (fname, e), outFile)

    quick = args.q
    bodies = args.p
    N = args.n
    runs = args.r
    benchmarks = args.b
    seq_dir = args.s
    mpi_dir = args.m
    processes = args.t

    tstart = datetime.datetime.now()
    run_benchmarks(quick, bodies, N, runs, benchmarks, seq_dir, mpi_dir, processes, outFile)
    delta = datetime.datetime.now() - tstart
    secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
    print("Total test time = %.2f seconds" % secs)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-q', metavar='quick', type=int, help='quick mode', default=0)
    parser.add_argument('-p', metavar='bodies', type=int, help='number of planet bodies', default=1000)
    parser.add_argument('-n', metavar='sim steps', type=int, help='number of simulation steps', default=200)
    parser.add_argument('-r', metavar='runs', type=int, help='number of runs', default=10)
    parser.add_argument('-f', metavar='outFile', type=str, help='output file name', default='YOIT2')
    parser.add_argument('-b', metavar='benchmarks', type=str, help='which benchmarks to run', default='123')
    parser.add_argument('-s', metavar='seq dir', type=str, help='dir of sequential executable', default='./gsim-seq/gsim')
    parser.add_argument('-bad', metavar='bad seq dir', type=str, help='dir of bad executable', default='./gsim-bad/gsim')
    parser.add_argument('-m', metavar='mpi dir', type=str, help='dir of mpi executable', default='./sim-mpi')
    parser.add_argument('-t', metavar='threads', type=int, help='process/thread count', default='1')
    args = parser.parse_args()

    run(args)

if __name__ == '__main__':
    main()
