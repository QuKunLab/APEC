#!/usr/bin/python
#
import numpy
import pandas
#import numba
from multiprocessing import Pool
import scipy.sparse
#
#
global GC_bias, peak_reads, nStep, sd
#
#
def mahalanobis_transform(aMatrix):
    (nSample, nFeature) = aMatrix.shape
    am_sum = aMatrix.sum(axis=0)
    ajaMatrix = numpy.zeros((nFeature, nFeature))
    for i in range(0, nFeature):
        for j in range(0, nFeature):
            ajaMatrix[i,j] = am_sum[i] * am_sum[j]
#    print ajaMatrix.sum(), aMatrix.sum()
    SaMatrix = numpy.dot(aMatrix.T, aMatrix) / float(nSample) - ajaMatrix / float(nSample**2)
#    print SaMatrix.sum(), nSample
    value, vector = numpy.linalg.eig(SaMatrix)
    value_inverseRoot = numpy.diag(-abs(value)**0.5)
    Sa_inverseRoot = numpy.dot(numpy.dot(vector, value_inverseRoot), vector.T)
    aMatrix_ave = aMatrix.mean(axis=0)
    aMatrix_ave = numpy.array([aMatrix_ave for i in range(0, nSample)])
    zMatrix = numpy.dot(Sa_inverseRoot, (aMatrix.T - aMatrix_ave.T))
    return zMatrix.T
#
#
#@numba.jit()
def single_sampling(par):
    GC_bias, peak_reads, nStep, sd, iIter = par[0], par[1], par[2], par[3], par[4]
    print 'permuted sampling ', iIter
    numpy.random.seed(12345+iIter)
    bias_step = (GC_bias.max() - GC_bias.min()) / float(nStep)
    read_step = (peak_reads.max() - peak_reads.min()) / float(nStep)
    sample = numpy.zeros(len(GC_bias),dtype=numpy.int)
    for ibias,bias_i in enumerate(GC_bias):
        bias_iIndex = int((bias_i - GC_bias.min()) // bias_step)
        read_iIndex = int((peak_reads[ibias] - peak_reads.min()) // read_step)
        bias_iIndex = min(nStep-1, max(0, bias_iIndex))
        read_iIndex = min(nStep-1, max(0, read_iIndex))
        peaks_inGrid = numpy.array([])
        ncount = 0
        while (len(peaks_inGrid)<=0) & (ncount<1000):
            ncount += 1
            bias_jIndex = bias_iIndex + int(numpy.rint(numpy.random.randn()*sd))
            read_jIndex = read_iIndex + int(numpy.rint(numpy.random.randn()*sd))
            while (bias_jIndex<0)|(bias_jIndex>=nStep):
                bias_jIndex = bias_iIndex + int(numpy.rint(numpy.random.randn()*sd))
            while (read_jIndex<0)|(read_jIndex>=nStep):
                read_jIndex = read_iIndex + int(numpy.rint(numpy.random.randn()*sd))
            bias_jStart = GC_bias.min() + bias_jIndex * bias_step
            read_jStart = peak_reads.min() + read_jIndex * read_step
            bias_jBin = numpy.where((bias_jStart<=GC_bias)&(GC_bias<bias_jStart+bias_step))[0]
            read_jBin = numpy.where((read_jStart<=peak_reads)&(peak_reads<read_jStart+read_step))[0]
            peaks_inGrid = numpy.intersect1d(bias_jBin, read_jBin)
        if ncount<1000:
            sample[ibias] = numpy.random.choice(peaks_inGrid)
        else:
            sample[ibias] = ibias
    print 'permuted sampling ', iIter, ' done!'
    return sample
#
#
def batch_sampling(GC_bias, peak_reads, nStep, sd, nIteration, np):
    matrix = numpy.vstack((GC_bias, peak_reads)).T
    matrix = mahalanobis_transform(matrix)
    GC_bias, peak_reads = matrix.T[0,:], matrix.T[1,:]
    kIterations = numpy.arange(0, nIteration, 1, dtype=int)
    parameters = []
    for iIter in kIterations:
        parameters.append([GC_bias, peak_reads, nStep, sd, iIter])
#
#    samples = []
#    for i in kIterations:
#        par = [GC_bias, peak_reads, nStep, sd, i]
#        sample = single_sampling(par)
#        samples.append(sample)
    pool = Pool(np)
    samples = pool.map(single_sampling, parameters)
    pool.close()
    pool.join()
    return numpy.array(samples)
#
#
#@numba.jit(nopython=True)
def expected_matrix(reads, TFmotif, GC_bias):
    (nCell, nPeak) = reads.shape   # X_matrix
    (nTF, nPeak) = TFmotif.shape   # M_matrix
    E_matrix = numpy.zeros(reads.shape)
    reads_sum = reads.sum()
    for iCell in range(0, nCell):
        reads_iCell = reads[iCell, :].sum()
        for jPeak in range(0, nPeak):
            E_matrix[iCell, jPeak] = reads[:, jPeak].sum() * reads_iCell / reads_sum
    return E_matrix
#
#
#
def deviation(MM, BB, XX, EE):
    MM = numpy.dot(MM, BB.T)
    YY = numpy.dot(MM, XX.T) - numpy.dot(MM, EE.T)
    denom = numpy.dot(MM, EE.T)
    YY = YY / denom
    return YY
#
#
def raw_deviation(TFmotif, reads, expected):
    (nCell, nPeak) = reads.shape
    (nTF, nPeak) = TFmotif.shape
    B_matrix = numpy.diag(numpy.ones(nPeak))
    raw_dev = deviation(TFmotif, B_matrix, reads, expected)
    return raw_dev
#
#
#
#
#
#
#
