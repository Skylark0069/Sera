import os

def main(data_save_dir, data_save_dir1):
    # 将 KU_net.xlsx 移动到缓存文件夹
    os.rename(data_save_dir + '\\'+ r'KU_net.xlsx', data_save_dir1 + '\\'+ r'KU_net.xlsx')
    # 将 真实结构差.xlsx 移动到缓存文件夹
    os.rename(data_save_dir + '\\'+ r'真实结构差.xlsx', data_save_dir1 + '\\'+ r'真实结构差.xlsx')
    # 将 结构差.xlsx 移动到缓存文件夹
    os.rename(data_save_dir + '\\'+ r'结构差.xlsx', data_save_dir1 + '\\'+ r'结构差.xlsx')

    # 将 参与成网边.xlsx 重命名为 Edge participating in structural identification.xlsx
    os.rename(data_save_dir + '\\'+ r'参与成网边.xlsx', data_save_dir + '\\'+ r'Edge participating in structural identification.xlsx')
    # 将 总表.xlsx 重命名为 SERA side table.xlsx
    os.rename(data_save_dir + '\\'+ r'总表.xlsx', data_save_dir + '\\'+ r'SERA side table.xlsx')
    # 将 网络传播.xlsx 重命名为 Total SERA network.xlsx
    os.rename(data_save_dir + '\\'+ r'网络传播.xlsx', data_save_dir + '\\'+ r'Total SERA network.xlsx')
    # 将 相似性.xlsx 重命名为 Similarity side table.xlsx
    os.rename(data_save_dir + '\\'+ r'相似性.xlsx', data_save_dir + '\\'+ r'Similarity side table.xlsx')

if __name__ == '__main__':
    main()