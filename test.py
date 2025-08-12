# import ENT3C --> runs ENT3C/__init__.py
# print(dir(ENT3C))

# import ENT3C.cli.main

# print(ENT3C.cli.main.__file__)
# print(ENT3C.cli.main.run_get_similarity)

# ENT3C_OUT = ENT3C.cli.main.run_get_entropy("config/test.sc.40kb.json")
# Similarity = ENT3C.cli.main.run_get_similarity("config/test.sc.40kb.json")
# if testing here, use pip -e install .
import ENT3C

print(ENT3C.__file__)

C_FN = "config/config.json"
group1 = "HFFc6"
group2 = "G401"

# ENT3C_OUT = ENT3C.run_get_entropy(C_FN)
# Similarity = ENT3C.run_get_similarity(C_FN)
ENT3C_OUT, Similarity = ENT3C.run_all(C_FN)


EUCLIDEAN = ENT3C.run_compare_groups(
    C_FN,
    group1,
    group2,
)

# print(EUCLIDEAN)
