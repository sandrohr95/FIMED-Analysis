class Config:
    __conf = {
        "GENE_EXPRESSION_DATA": "data_dir/Gene_expression_data_dir/",
        "GRN_PANDA_DIR": '/home/antonio/Anaconda_Projects/Python_Projects/FIMED-Analysis/Gene_Regulation_Networks/data_dir/Results/',
        "GRN_PANDA_ENV": '/home/antonio/Anaconda_Projects/Python_Projects/FIMED-Analysis/pypanda/pypandaenv/bin/python',
        "GRN_PANDA_SCRIPT": '/home/antonio/Anaconda_Projects/Python_Projects/FIMED-Analysis/pypanda/run_panda.py'

    }

    @staticmethod
    def config(name):
        return Config.__conf[name]

