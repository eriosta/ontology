from owlready2 import get_ontology, ThingClass
from collections import defaultdict
import json

ONTOLOGY_URL = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/doid.owl"
CANCER_ROOT = "DOID:162"
CELLULAR_PROLIF_ROOT = "DOID:14566"
ORGAN_SYSTEM_CANCER = "DOID:0050686"

print("🔄 Loading ontology...")
onto = get_ontology(ONTOLOGY_URL).load()

label_map = {}
parent_map = defaultdict(list)
child_map = defaultdict(list)

print("🔍 Indexing ontology...")
for cls in onto.classes():
    if "DOID_" not in cls.name:
        continue
    curie = "DOID:" + cls.name.split("_")[-1]
    label_map[curie] = cls.label.first() or cls.name

    for parent in cls.is_a:
        if isinstance(parent, ThingClass) and hasattr(parent, "name") and "DOID_" in parent.name:
            parent_curie = "DOID:" + parent.name.split("_")[-1]
            parent_map[curie].append(parent_curie)
            child_map[parent_curie].append(curie)

# ✅ Step: Recursive path tracing
def trace_paths_to_root(curie, path=None):
    path = [curie] + (path or [])
    parents = parent_map.get(curie, [])
    if not parents:
        return [path]
    paths = []
    for p in parents:
        paths.extend(trace_paths_to_root(p, path))
    return paths

# ✅ Step: Find leaf nodes (no children)
print("🌿 Finding leaf nodes in cancer subtree...")
leaf_nodes = []
for curie in label_map:
    # is leaf and in cancer lineage
    if curie not in child_map:
        paths = trace_paths_to_root(curie)
        if any(CANCER_ROOT in p and CELLULAR_PROLIF_ROOT in p and ORGAN_SYSTEM_CANCER in p for p in paths):
            leaf_nodes.append(curie)

print(f"🌿 Found {len(leaf_nodes)} cancer leaf terms")

# ✅ Step: Construct dictionary from leaves
leaf_hierarchy = {}
for curie in leaf_nodes:
    paths = trace_paths_to_root(curie)
    paths = [p for p in paths if CANCER_ROOT in p and CELLULAR_PROLIF_ROOT in p and ORGAN_SYSTEM_CANCER in p]
    if not paths:
        continue
    label_paths = [[label_map.get(c, c) for c in path] for path in paths]
    leaf_hierarchy[curie] = {
        "label": label_map[curie],
        "paths_to_root": paths,
        "label_paths_to_root": label_paths
    }

# ✅ Step: Save result
import os
os.makedirs("dictionaries/disease", exist_ok=True)
output_file = "dictionaries/disease/doid_cancer_leaf_paths.json"
with open(output_file, "w") as f:
    json.dump(leaf_hierarchy, f, indent=2)

print(f"✅ Saved {len(leaf_hierarchy)} leaf nodes to {output_file}")