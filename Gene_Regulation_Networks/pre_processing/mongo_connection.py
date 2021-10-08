#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 17:03:43 2018

@author: Sandro Hurtado
"""

from pymongo import MongoClient
from bson.objectid import ObjectId
import gridfs


class MongoConnection:

    def __init__(self):
        mongo_client = MongoClient("mongodb://root:kha0sd3v@192.168.213.38:27017")
        # mongo_client = MongoClient("mongodb://root:root@localhost:27017")
        self.db = mongo_client.FIMED
    
    def read_from_mongo(self, input):
        fs = gridfs.GridFS(self.db, collection='Clinical_Samples')
        
        for grid_out in fs.find({"_id": ObjectId(input)},no_cursor_timeout=True):
            data = grid_out.read().decode("utf-8")
            return data
        
    def get_sample_name(self, input):
        query = self.db.Clinical_Samples.files.find_one({"_id": ObjectId(input)})
        return query['filename']
    
    
    
        
    

    
    
    
    

   
    



