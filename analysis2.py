import numpy as np
import matplotlib.pyplot as plt
import csv

from config import simInputData

class simAnalysisData:
    pressure = []
    oxygen = []
    vegf = []
    signal = []

def collect_data(sad:simAnalysisData, sid:simInputData, in_nodes, pnow, vnow, snow):
    data = [pnow[in_nodes[0]], np.average(vnow), np.average(snow)]
    sad.pressure.append(data[0])
    sad.vegf.append(data[1])
    sad.signal.append(data[2])
    f = open(sid.dirname+'/data.csv', 'a')
    f.write(f'{data[0]} {data[1]} {data[2]} \r')
    f.close()

def plot_data(sid:simInputData):
    f = open(sid.dirname+'/data.csv', 'r')
    data = csv.reader(f, delimiter = ' ')
    pressure = []
    vegf = []
    signal = []
    for row in data:
        pressure.append(row[0])
        vegf.append(row[1])
        signal.append(row[2])
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(pressure)
    ax2.plot(vegf)
    ax3.plot(signal)
    plt.savefig(sid.dirname + "/data.png")
    plt.close()