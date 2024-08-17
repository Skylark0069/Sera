# 将csv文件转成xlsx文件
import csv
import os
import sys
from openpyxl import Workbook
def csv2xlsx(csv_file, xlsx_file):
    wb = Workbook()
    ws = wb.active
    with open(csv_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            ws.append(row)
    # 表格的第一列和第二列数据类型转化为int, 第三列转化为float, 第四列转化为str
    for row in ws.iter_rows(min_row=1, max_col=4):
        row[0].value = int(row[0].value)
        row[1].value = int(row[1].value)
        row[2].value = float(row[2].value)
        row[3].value = str(row[3].value)
    wb.save(xlsx_file)
    # 删除csv文件
    os.remove(csv_file)
        
def main(data_save_dir):
    csv2xlsx(data_save_dir + '\\'+ r'结构差.csv', data_save_dir + '\\'+ r'结构差.xlsx')

if __name__ == '__main__':
    main()
    