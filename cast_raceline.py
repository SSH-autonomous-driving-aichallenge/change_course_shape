import csv
import os
import math

def convert_vehicle_data_csv(input_filepath, output_filepath,
                             a_lat_max=14.0,  # 許容横加速度 [m/s²]
                             v_max=60.0       # 最高速度 [m/s]
                            ):
    """
    車両データを指定された形式のCSVに変換します。
    速度は曲率に対応して計算します。

    Args:
        input_filepath (str): 入力CSVファイルのパス。
        output_filepath (str): 出力CSVファイルのパス。
        a_lat_max (float): 許容横加速度 [m/s²]。
        v_max (float): 曲率ゼロ時の最高速度 [m/s]。
    """
    output_header = ['x', 'y', 'z', 'x_quat', 'y_quat', 'z_quat', 'w_quat', 'speed']

    with open(input_filepath, 'r', newline='', encoding='utf-8') as infile, \
         open(output_filepath, 'w', newline='', encoding='utf-8') as outfile:

        reader = csv.reader(infile, delimiter=',')
        writer = csv.writer(outfile, delimiter=',')
        writer.writerow(output_header)
        next(reader, None)  # ヘッダー行をスキップ

        for row in reader:
            if not row:
                continue

            try:
                x = row[0].strip()
                y = row[1].strip()
                kappa = float(row[2].strip())

                # 曲率から速度を計算（横加速度制限式）
                if abs(kappa) > 1e-4:
                    v = math.sqrt(a_lat_max / abs(kappa))
                else:
                    v = v_max  # ほぼ直線なら最高速度
                # 最高速度をオーバーしないようにクリップ
                speed = min(v, v_max)

                new_row = [
                    x, y,
                    0,  # z
                    0, 0, 0, 0,  # quaternion（固定）
                    f"{speed:.3f}"
                ]
                writer.writerow(new_row)

            except (IndexError, ValueError):
                # データ不備行はスキップ
                continue

    print(f"✅ 変換が完了しました: {output_filepath}")

# 実行部はそのまま
if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(script_dir, 'traj.csv')
    output_path = os.path.join(script_dir, 'raceline_awsim_15km.csv')
    convert_vehicle_data_csv(input_path, output_path)

