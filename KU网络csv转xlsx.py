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
    for row in ws.iter_rows(min_row=1, max_col=6):
        if row[0].value != 'nist':
            row[0].value = int(row[0].value)
        else:
            row[0].value = str(row[0].value)
        row[1].value = str(row[1].value)
        if row[2].value != 'nist':
            row[2].value = int(row[2].value)
        else:
            row[2].value = str(row[2].value)
        row[3].value = str(row[3].value)
        row[4].value = float(row[4].value)
        row[5].value = str(row[5].value)
    wb.save(xlsx_file)
    # 删除csv文件
    os.remove(csv_file)
        
def main(data_save_dir):
    csv2xlsx(data_save_dir + '\\'+ r'KU_net.csv', data_save_dir + '\\'+ r'KU_net.xlsx')

if __name__ == '__main__':
    main()

    