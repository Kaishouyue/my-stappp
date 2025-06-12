import re
import numpy as np
import os

debug = 0  # Debug switch

# Get user input for filename (without path and extension)
user_file = input("Please enter the filename (without path and extension, e.g. test):")
filename = fr'c:\Users\lee\Desktop\STAPpp-master\src\data\{user_file}.out'


# 1. Read element type
with open(filename, 'r') as f:
    lines = f.readlines()

# --- Section location ---
node_start = elem_start = disp_start = None
for i, line in enumerate(lines):
    if 'N O D A L   P O I N T   D A T A' in line:
        node_start = i+3  # Skip 3 header lines
    if 'E L E M E N T   I N F O R M A T I O N' in line:
        elem_start = i+2  # Skip 2 header lines
    if 'D I S P L A C E M E N T S' in line:
        disp_start = i+2  # Skip 2 header lines
    if 'ELEMENT TYPE' in line:
        m = re.search(r'=\s*(\d+)', line)
        if m:
            element_type = int(m.group(1))

if debug == 1:
    print(f"element_type = {element_type}")

# --- Switch-case style for reading node coordinates ---
def read_nodes_q4(lines, node_start):
    node_coords = []
    for line in lines[node_start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+', line):
            parts = line.split()
            x, y = map(float, parts[-2:])
            node_coords.append([x, y])
        elif line.strip() == '':
            break
    return np.array(node_coords)

def read_nodes_truss(lines, node_start):
    node_coords = []
    for line in lines[node_start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+', line):
            parts = line.split()
            x, y, z = map(float, parts[-3:])
            node_coords.append([x, y, z])
        elif line.strip() == '':
            break
    return np.array(node_coords)

read_nodes_switch = {
    2: read_nodes_q4,
    1: read_nodes_truss,
}

node_coords = read_nodes_switch.get(element_type, lambda l, s: np.array([]))(lines, node_start)

if debug == 1:
    print("node_coords =")
    print(node_coords)

# --- Switch-case style for reading element connectivity ---
def read_elem_connect_q4(lines, elem_start):
    elem_connect = []
    for line in lines[elem_start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+', line):
            parts = line.split()
            elem_connect.append([int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])])
        elif line.strip() == '':
            break
    return elem_connect

def read_elem_connect_truss(lines, elem_start):
    elem_connect = []
    for line in lines[elem_start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+\s+\d+', line):
            parts = line.split()
            elem_connect.append([int(parts[1]), int(parts[2])])
        elif line.strip() == '':
            break
    return elem_connect

read_elem_connect_switch = {
    2: read_elem_connect_q4,
    1: read_elem_connect_truss,
}

elem_connect = read_elem_connect_switch.get(element_type, lambda l, s: [])(lines, elem_start)

if debug == 1:
    print("elem_connect =")
    print(elem_connect)

# --- Switch-case style for reading displacements ---
def read_disp_q4(lines, disp_start):
    displacements = []
    for line in lines[disp_start:]:
        if re.match(r'\s*\d+\s+[-\deE.+]+\s+[-\deE.+]+\s+[-\deE.+]+', line):
            parts = line.split()
            dx = float(parts[1])
            dy = float(parts[2])
            displacements.append([dx, dy])
        elif line.strip() == '':
            break
    return np.array(displacements)

def read_disp_truss(lines, disp_start):
    displacements = []
    for line in lines[disp_start:]:
        if re.match(r'\s*\d+\s+[-\deE.+]+\s+[-\deE.+]+\s+[-\deE.+]+', line):
            parts = line.split()
            dx = float(parts[1])
            dy = float(parts[2])
            dz = float(parts[3]) if len(parts) > 3 else 0.0
            displacements.append([dx, dy, dz])
        elif line.strip() == '':
            break
    return np.array(displacements)

read_disp_switch = {
    2: read_disp_q4,
    1: read_disp_truss,
}

displacements = read_disp_switch.get(element_type, lambda l, s: np.array([]))(lines, disp_start)

if debug == 1:
    print("displacements =")
    print(displacements)

# 计算L2范数
if element_type == 2:
    # Q4单元，二维节点
    u_exact = 0.12 * node_coords[:,0] * (node_coords[:,1] - 0.5)
    v_exact = -0.06 * node_coords[:,0]**2
    u_num = displacements[:,0]
    v_num = displacements[:,1]
    for i in range(len(node_coords)):
        print(f"节点{i+1}: 计算位移(u, v)=({u_num[i]:.6e}, {v_num[i]:.6e}), 精确位移(u, v)=({u_exact[i]:.6e}, {v_exact[i]:.6e})")
    error = np.sqrt(np.sum((u_num - u_exact)**2 + (v_num - v_exact)**2) / len(node_coords))
    print(f"L2范数（均方根误差）: {error}")