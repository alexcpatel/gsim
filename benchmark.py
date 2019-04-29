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

def doRun(cmd_args, benchmark, quick, outFile):
    cmdline = " ".join(cmd_args)
    tstart = datetime.datetime.now()

    try:
        print("Running '%s'" % cmdline)
        simProcess = subprocess.Popen(cmd_args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        stream = simProcess.communicate()
        if quick:
            stream = None
        else:
            stream = stream[0].splitlines()
        simProcess.wait()
        returnCode = simProcess.returncode
    except Exception as e:
        print("Execution of command '%s' failed with exception: %s" % (cmdline, e))
        return None, None
    if returnCode == 0:
        delta = datetime.datetime.now() - tstart
        secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
        return secs, stream
    else:
        print("Execution of command '%s' failed with code: %s" % (cmdline, returnCode))
        return None, None

def bestRun(cmd_args, benchmark, runs, quick, outFile):
    sofar = sys.maxsize
    best_stream = None
    for r in range(runs):
        if runs > 1:
            print("Run #%d" % (r+1))
        secs, stream = doRun(cmd_args, benchmark, quick, outFile)
        if quick:
            assert(stream == None)
        if secs == None:
            return None
        if secs < sofar:
            sofar, best_stream = secs, stream
    print(sofar, best_stream)
    return sofar, best_stream

def run_benchmark(quick, clusters, bodies, N, runs, benchmark, seq_dir, mpi_dir, processes, outFile, theta):
    cmd_args = [mpi_dir, str(clusters), str(bodies), str(N), str(theta), str(seed)]
    secs, stream = bestRun(cmd_args, benchmark, runs, quick, outFile)
    if quick:
        assert(stream == None)

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
        if bodies < 4:
            clusters = 1
            return clusters, bodies
        else:
            clusters = math.sqrt(bodies)
            return clusters, bodies

def streams_equal(s1, s2, e):
    for position in range(len(s1)):
        p1 = s1[position]
        p2 = s2[position]
        p1 = p1.split()
        p2 = p2.split()
        if abs(p1[-1] - p2[-1]) > e:
            return False
        if abs(p1[-2] - p2[-2]) > e:
            return False
    return True

def run_benchmarks(quick, bodies, N, runs, benchmarks, seq_dir, mpi_dir, processes, thetas, epsilon, outFile):
    if processes == 1:
        quick = True
    for process in range(processes):
        for theta in thetas.split(','):
            outmsg("PARAMETERS ----- bodies: %d, steps: %d, runs: %d, processes: %d, theta: %f" %
                    (bodies, N, runs, process, float(theta)), outFile)
            for benchmark in list(benchmarks):
                clusters, bodies = get_benchmark(benchmark, bodies)
                if processes > 1:
                    results_mpi = run_benchmark(quick, clusters, bodies, N, runs, benchmark, seq_dir, mpi_dir,
                                            processes, outFile, float(theta))
                    outmsg("Benchmark MPI %s -- Seconds: %f, NPM: %f" % (benchmark, results_mpi[0], results_mpi[1]), outFile)
                    if not streams_equal(results_mpi[2], results_seq[2], epsilon):
                        outmsg("Implementations did NOT produce equivalent results!!")
                        return -1

def run(args):
    fname = 'benchmark_out/' + args.f
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
    epsilon = args.e
    thetas = args.thetas

    tstart = datetime.datetime.now()
    run_benchmarks(quick, bodies, N, runs, benchmarks, seq_dir, mpi_dir, processes, thetas, epsilon, outFile)
    delta = datetime.datetime.now() - tstart
    secs = delta.seconds + 24 * 3600 * delta.days + 1e-6 * delta.microseconds
    print("Total test time = %.2f seconds" % secs)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-q', metavar='quick', type=int, help='quick mode', default=1)
    parser.add_argument('-p', metavar='bodies', type=int, help='number of planet bodies', default=1000)
    parser.add_argument('-n', metavar='sim steps', type=int, help='number of simulation steps', default=200)
    parser.add_argument('-r', metavar='runs', type=int, help='number of runs', default=10)
    parser.add_argument('-f', metavar='outFile', type=str, help='output file name', default='temp_log')
    parser.add_argument('-b', metavar='benchmarks', type=str, help='which benchmarks to run', default='123')
    parser.add_argument('-s', metavar='seq dir', type=str, help='dir of sequential executable', default='./gsim-seq/gsim-seq')
    parser.add_argument('-m', metavar='mpi dir', type=str, help='dir of mpi executable', default='./gsim-omp/')
    parser.add_argument('-t', metavar='threads', type=int, help='process/thread count', default='8')
    parser.add_argument('-e', metavar='epsilon', type=float, help='epsilon wiggle room when testing accuracy', default=1000.0)
    parser.add_argument('-thetas', metavar='thetas', type=str, help='list of the thetas to be used for the benchmarking', default='0.0, 0.5, 1.0')
    args = parser.parse_args()

    run(args)

if __name__ == '__main__':
    main()
