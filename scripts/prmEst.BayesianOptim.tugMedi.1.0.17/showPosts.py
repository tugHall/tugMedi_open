#!/bin/python3

# coding: shift-jis
import sys
import logging
import argparse


#logging.basicConfig(level=logging.DEBUG)


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def get_top_replicates(nn, file):

    with open(file, "r", encoding="us-ascii") as ff:
        row = ff.read()
    
    lines = row.strip().split("\n")
    logging.debug(lines)

    data        = [(line.split("\t")[0], float(line.split("\t")[1])) for line in lines]
    sorted_data = sorted(data, key=lambda x: x[1])  
    res = {sorted_data[ii][0]: sorted_data[ii][1] for ii in range(min(nn, len(sorted_data)))} 

    return res
    
    
def process_file(prm1_path, prm2_path):

    with open(prm1_path, "r", encoding="us-ascii") as ff:
        lines = ff.readlines()
    prms = {
        line.split("\t")[0]: line.split("\t")[1].strip()
        for line in lines if line.strip() and "\t" in line
    }
    
    with open(prm2_path, "r", encoding="us-ascii") as ff:
        lines = ff.readlines()
    EFtimes = {
        line.split("\t")[1]: line.split("\t")[0]
        for line in lines[1:] if line.strip() and "\t" in line
    }
    
    ret = {"m0": prms.get("m0", None), # m0 and dN among all
           "dN": prms.get("dN", None)}
    ret.update(EFtimes)
    
    return ret


def get_toDisplay(mm, dir_baseReps, prm1, prm2):
    
    results = []
    for num, dist in mm:
        
        dir_path  = f"{dir_baseReps}/{str(num).zfill(6)}" # 0-padding, 6 digits 
        prm1_path = f"{dir_path}/{prm1}"
        prm2_path = f"{dir_path}/{prm2}"
        
        ret = process_file(prm1_path, prm2_path)
        ret.update({ "rep": num, "ABCdist": dist })
        results.append(ret)
    return results


def display(results):

    ordered_keys = ["rep", "ABCdist", "m0", "dN"] # keys to show
    remaining_keys = sorted(
        {key for ret in results for key in ret.keys() if key not in ordered_keys}
    )
    final_order = ordered_keys + remaining_keys

    header = "\t".join(final_order)
    print(header)
    
    for ret in results:
        row = "\t".join(
            [f"{ret.get(key, ''):.5g}"  # 5 digits float
                if key == "ABCdist" and isinstance( ret.get(key), (int, float) ) 
                else str( ret.get(key, "") )
             for key in final_order]
        )
        print(row)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


prm1 = "Input/parameters.txt"
prm2 = "Input/EF.Rint.txt"


def main():

# nn = 10
# file_ABCdist = "./ABC2/TCGA-55-7903-01A-11D-2167-08/ABC_dist_tumor_content=0.80.txt"
# dir_baseReps = "./work/TCGA-55-7903-01A-11D-2167-08"

    parser = argparse.ArgumentParser(description="Get posteriors selected by ABC.")
    
    parser.add_argument("--nn",           type=int, required=True, help="# of top replicates with shortest distance (eg, 10)")
    parser.add_argument("--file_ABCdist", type=str, required=True, help="Path to the ABC distance file")
    parser.add_argument("--dir_baseReps", type=str, required=True, help="Path to base directory to replicates")

    args = parser.parse_args()
    
    nn = args.nn
    file_ABCdist = args.file_ABCdist
    dir_baseReps = args.dir_baseReps

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    val = get_top_replicates(nn, file_ABCdist)

    mm = sorted(val.items(), key=lambda x: x[1])
    logging.debug("mm:", mm)

    ret = get_toDisplay(mm, dir_baseReps, prm1, prm2)
    display(ret)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if __name__ == "__main__":
    main()

