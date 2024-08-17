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

    for row in ws.iter_rows(min_row=2, min_col=1, max_col=7):
        row[0].value, row[2].value, row[6].value = int(row[0].value), int(row[2].value), int(row[6].value)
        row[4].value = float(row[4].value)
        row[1].value, row[3].value, row[5].value = str(row[1].value), str(row[3].value), str(row[5].value)

    wb.save(xlsx_file)
    # 删除csv文件
    os.remove(csv_file)
        
def main(data_save_dir, data_save_name):
    csv2xlsx(data_save_dir + '\\'+ r'网络传播.csv', data_save_dir + '\\'+ r'网络传播.xlsx')

if __name__ == '__main__':
    main()