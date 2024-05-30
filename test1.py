import math
import matplotlib.pyplot as plt


# 定义一个空列表来存储原子数据
atoms = [[], [], []]
file_names = ["water_273.trj", "water_300.trj", "water_400.trj"]

# 打开数据文件
for file_index, file_name in enumerate(file_names):
    with open(file_name, "r") as file:
        for line in file:
            # 去掉行尾的换行符
            line = line.strip()
            if "ITEM: TIMESTEP" in line:
                line = file.readline().strip()
                timestep = int(int(line) / 1000)
                # print(timestep)
                atoms[file_index].append([])
                for i in range(7):
                    file.readline()
                # 读取一帧原子数据
                while True:
                    atom_line = file.readline().strip()
                    # print(atom_line)
                    atom_data = atom_line.split()
                    atom = {
                        "id": int(atom_data[0]),
                        "mol": int(atom_data[1]),
                        "xu": float(atom_data[2]),
                        "yu": float(atom_data[3]),
                        "zu": float(atom_data[4]),
                        "type": int(atom_data[5]),
                        "element": atom_data[6],
                    }
                    atoms[file_index][timestep].append(atom)
                    if atom["id"] == 4500:
                        break
                # 一旦读取完所有原子数据，就退出循环

# 打印原子数据
print(atoms[0][0][0])
# print(atoms[1][4513509-11]import math

# atoms为指定温度
atoms_273 = atoms[0]
atoms_300 = atoms[1]
atoms_400 = atoms[2]


def get_density_distribution(atoms_t, dr, sframe, eframe, mtype, stype):
    density_frame = []
    density = []
    mnum = 0
    bintotal = 100000000
    max_distance = 0
    for i in range(0, bintotal + 1):
        density.append(0)
        density_frame.append(0)
    for i in range(sframe, eframe + 1):
        for j in range(0, 4500):
            mnum += 1
            if atoms_t[i][j]["type"] == mtype:
                for k in range(0, 4500):
                    if atoms_t[i][k]["type"] == stype:
                        distance = math.sqrt(
                            (atoms_t[i][j]["xu"] - atoms_t[i][k]["xu"]) ** 2
                            + (atoms_t[i][j]["yu"] - atoms_t[i][k]["yu"]) ** 2
                            + (atoms_t[i][j]["zu"] - atoms_t[i][k]["zu"]) ** 2
                        )
                        if distance > max_distance:
                            max_distance = distance
                        nbin = int(distance / dr)
                        density_frame[nbin] += 1
                        #print("nbin", nbin)
                        #print("density_frame[nbin]", density_frame[nbin])
        bintotal = int(max_distance / dr)
        density_frame = [x / mnum for x in density_frame]
        print(density_frame[:bintotal])
        for bin in range(0, bintotal):
            # print(density_frame[bin])
            density_frame[bin] = density_frame[bin] / (
                (4 / 3) * math.pi * dr * ((bin + 1) ** 3 - bin**3)
            )
            density[bin] += density_frame[bin]
            # print(density_frame[bin])
            # print('###')
        density_frame[bintotal] = density_frame[bintotal] / (
            (4 / 3) * math.pi * (max_distance**3 - (dr * bintotal) ** 3)
        )
        # print(bintotal)
        density[bintotal] += density_frame[bintotal]
        # density=density[:bintotal]
        print("####")
        print(density[:bintotal])
    density = [x / (eframe - sframe + 1) for x in density]
    # print(bintotal)
    return density, max_distance


def pic_line_chart(L,max_dist,dr):
    nbin=int(max_dist/dr)
    x=[i*dr for i in range(0,nbin+1)]
    y=L[:nbin+1]
    plt.plot(x,y)
    plt.savefig("density_distribution.png")

#task1求解
dr=1e-1
density,max_dist=get_density_distribution(atoms_273,dr,0,0,1,1)
density=list(density)
pic_line_chart(density,max_dist,dr)
    