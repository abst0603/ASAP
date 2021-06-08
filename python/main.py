from extract_cycles import extract_cycles
from Voting import voting
import pickle

# extracted_cycles = extract_cycles('ASAP100/samples180_6e-6_0.20.mat', 3, 0.1)
file_name = "extracted_cycles.pkl"

# open_file = open(file_name, "wb")
# pickle.dump(extracted_cycles, open_file)
# open_file.close()

open_file = open(file_name, "rb")
loaded_list = pickle.load(open_file)
open_file.close()

candidates_list = voting(loaded_list,2,0.4)# radius is chosen twice the radius of sampling

print('done')
