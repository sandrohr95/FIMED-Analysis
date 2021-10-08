import os
import pre_processing.parse as pr
import numpy as np
import pre_processing.normalize as nor
import pre_processing.mongo_connection as mongo

# path = '/home/antonio/Anaconda_Projects/Python_Projects/Gene_Regulation_Networks'
# os.chdir(path)


def run_norm_mongodb(id_list):
    datalist = []
    flags = []
    genes = []

    mongodb = mongo.MongoConnection()
    for id in id_list:
        dt = mongodb.read_from_mongo(id)
        dt = pr.parse(dt)
        if (len(genes) > 0):
            if (np.all(np.array(dt[0]) == np.array(genes)) != True):
                print('Gene lists do not match')
                exit()
        genes = dt[0]  # Gene name
        datalist.append(dt[2])  # Gene values
        flags.append(dt[3])

    return nor.process(genes, dt[1], id_list, datalist, flags)
