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

def doRun(cmd_args, benchmark, outFile, quick):
    cmdline = " ".join(cmd_args)
    tstart = datetime.datetime.now()
    stream = []
    try:
        print("Running '%s'" % cmdline)
        simProcess = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        simProcess.wait()
        returnCode = simProcess.returncode
        if not quick:
            stream = simProcess.communicate()
    except Exception as e:
        print("Execution of command '%s' failed with exception: %s" % (cmdline, e))
        return None
    if returnCode == 0:
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        return secs, stream
    else:
        print("Execution of command '%s' failed with code: %s" % (cmdline, returnCode))
        return None

def bestRun(cmd_args, benchmark, runs, outFile, quick):
    sofar = [sys.maxsize]
    sostream = []
    for r in range(runs):
        if runs > 1:
            print("Run #%d" % (r+1))
        secs, stream = doRun(cmd_args, benchmark, outFile, quick)
        if secs == None:
            return None
        sofar, sostream = sofar, sostream if sofar < secs else secs, stream
    return sofar, sostream

def run_benchmark(quick, clusters, bodies, N, runs, benchmark, seq_dir, mpi_dir, processes, outFile):
    if processes > 1:
        cmd_args = [] # TODO
    else:
        cmd_args = [seq_dir, str(clusters), str(bodies), str(N), str(seed)]
    secs, stream = bestRun(cmd_args, benchmark, runs, outFile, quick)

    if not secs == None:
        bmoves = bodies * N
        npm = 1e9 * secs/bmoves
        return secs, npm, stream
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
        clusters = bodies // 4 if bodies > 8 else 2
        return clusters, bodies

def stream_equal(stream1, stream2):
    stream1 = stream1.splitlines()[1:]
    stream2 = stream2.splitlines()[1:]
    if len(stream1) != len(stream2):
        return False
    for i in range(len(stream1)):
        if stream1[i] != stream2[i]:
            return False
    return True

def run_benchmarks(quick, bodies, N, runs, benchmarks, seq_dir, mpi_dir, processes, outFile):
    tcount = 0
    rcount = 0
    sum = 0.0
    refSum = 0.0

    if processes == 1:
        outmsg("PARAMETERS ----- bodies: %d, steps: %d, runs: %d, processes: %d" %
                (bodies, N, runs, processes), outFile)
        for benchmark in list(benchmarks):
            clusters, bodies = get_benchmark(benchmark, bodies)
            results = run_benchmark(quick, clusters, bodies, N, runs, benchmark, seq_dir, mpi_dir, processes, outFile)
            outmsg("Benchmark %s -- Seconds: %f, NPM: %f" % (benchmark, results[0], results[1]), outFile)
    else: #we're parallel bois
        for benchmark in list(benchmarks):
            clusters, bodies = get_benchmarks(benchmark, bodies)
            results_seq = run_benchmark(quick, clusters, bodies, N,
                            runs, benchmark, seq_dir, mpi_dir, processes=1, outFile)
            outmsg("Benchmark %s -- Seconds: %f, NPM: %f" % (benchmark, results[0], results[1]), outFile)
            results_par = run_benchmark(quick, clusters, bodies, N,
                            runs, benchmark, seq_dir, mpi_dir, processes, outFile)
            outmsg("Benchmark %s -- Seconds: %f, NPM: %f" % (benchmark, results[0], results[1]), outFile)
            if not quick: # do error checking
                if not streams_equal(results_seq[2], results_par[2]):
                    print("ERROR: Parallel Implementation is NOT CORRECT")
                    outmsg("PARALLEL IMPLEMENTATION NOT CORRECT")

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
    parser.add_argument('-f', metavar='outFile', type=str, help='output file name', default='YOIT')
    parser.add_argument('-b', metavar='benchmarks', type=str, help='which benchmarks to run', default='123')
    parser.add_argument('-s', metavar='seq dir', type=str, help='dir of sequential executable', default='./gsim')
    parser.add_argument('-m', metavar='mpi dir', type=str, help='dir of mpi executable', default='./sim-mpi')
    parser.add_argument('-t', metavar='threads', type=int, help='process/thread count', default='1')
    args = parser.parse_args()

    run(args)

if __name__ == '__main__':
    main()
