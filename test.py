# import ENT3C --> runs ENT3C/__init__.py
# print(dir(ENT3C))

# import ENT3C.cli.main

# print(ENT3C.cli.main.__file__)
# print(ENT3C.cli.main.run_get_similarity)

# ENT3C_OUT = ENT3C.cli.main.run_get_entropy("config/test.sc.40kb.json")
# Similarity = ENT3C.cli.main.run_get_similarity("config/test.sc.40kb.json")
import ENT3C

print(ENT3C.__file__)

print("hello")
# ENT3C_OUT = ENT3C.run_get_entropy("config/test.sc.40kb.json")
# Similarity = ENT3C.run_get_similarity("config/test.sc.40kb.json")
ENT3C_OUT, Similarity = ENT3C.run_all("config/config.test.json")
print(ENT3C_OUT)
print(Similarity)
