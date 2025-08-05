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

# ENT3C_OUT = ENT3C.run_get_entropy("config/test.sc.40kb.json")
# Similarity = ENT3C.run_get_similarity("config/test.sc.40kb.json")

# ENT3C_OUT, Similarity = ENT3C.run_all("config/.json")


group1 = "ODC_AD"
group2 = "ODC_CTL"

Similarity = ENT3C.run_get_similarity("./config/config.10kb.json")
DIFFERENCE = ENT3C.run_locate_largest_euclidean_diff(
    "./config/AD.CTL.config.10kb.json", group1, group2
)

Similarity = ENT3C.run_get_similarity("./config/config.5kb.json")
DIFFERENCE = ENT3C.run_locate_largest_euclidean_diff(
    "./config/AD.CTL.config.5kb.json", group1, group2
)
