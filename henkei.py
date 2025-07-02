import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import messagebox
import numpy as np

class WaypointEditor:
    """
    CSVファイルから読み込んだウェイポイントを編集するGUIアプリケーション。

    主な機能:
    - CSVファイルを読み込み、ウェイポイントを可視化する。
    - 2点をクリックして、その間のウェイポイント（周回データとして最短経路）を一括選択する。
    - 選択した点をターゲット（右クリック地点）に向けて、指定した割合で移動させる。
    - 選択した点とその周辺の点に対して平滑化（移動平均）を適用する。
    - 変更をCSVファイルに保存する。
    """
    def __init__(self, master, filepath):
        self.master = master
        self.filepath = filepath
        self.selected_indices = []
        self.selection_endpoints = [] # 選択の始点と終点を格納

        # --- データ読み込み ---
        try:
            self.df = pd.read_csv(self.filepath)
            if not {'x', 'y'}.issubset(self.df.columns):
                raise ValueError("CSV must contain 'x' and 'y' columns.")
        except FileNotFoundError:
            messagebox.showerror("エラー", f"ファイルが見つかりません: {self.filepath}")
            self.master.quit()
            return
        except Exception as e:
            messagebox.showerror("エラー", f"ファイルの読み込みに失敗しました: {e}")
            self.master.quit()
            return

        # --- GUIのセットアップ ---
        self.master.title("ウェイポイントエディタ")
        self.master.geometry("1000x800")

        # --- トップフレーム（操作パネル） ---
        top_frame = tk.Frame(self.master)
        top_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)

        self.save_button = tk.Button(top_frame, text="CSVに保存", command=self.save_csv)
        self.save_button.pack(side=tk.LEFT, padx=5)

        tk.Label(top_frame, text="移動率:").pack(side=tk.LEFT, padx=(10, 0))
        self.move_factor = tk.DoubleVar(value=0.3)
        self.slider = tk.Scale(top_frame, from_=0.0, to=1.0, resolution=0.05,
                               orient=tk.HORIZONTAL, variable=self.move_factor)
        self.slider.pack(side=tk.LEFT, padx=5)

        # --- 平滑化コントロール ---
        smoothing_frame = tk.LabelFrame(top_frame, text="平滑化")
        smoothing_frame.pack(side=tk.LEFT, padx=10)

        tk.Label(smoothing_frame, text="窓サイズ(±):").pack(side=tk.LEFT, padx=(5,0))
        self.smooth_window = tk.IntVar(value=2)
        tk.Spinbox(smoothing_frame, from_=1, to=10, width=5, textvariable=self.smooth_window).pack(side=tk.LEFT, padx=(0,5))

        tk.Label(smoothing_frame, text="周辺点:").pack(side=tk.LEFT, padx=(5,0))
        self.smooth_padding = tk.IntVar(value=1)
        tk.Spinbox(smoothing_frame, from_=0, to=10, width=5, textvariable=self.smooth_padding).pack(side=tk.LEFT, padx=(0,5))
        
        self.smooth_button = tk.Button(smoothing_frame, text="平滑化を適用", command=self.smooth_selection)
        self.smooth_button.pack(side=tk.LEFT, padx=5, pady=5)

        instructions = "点2つを左クリックで区間選択。右クリックで選択点をカーソル方向へ移動。"
        self.label = tk.Label(top_frame, text=instructions, fg="blue")
        self.label.pack(side=tk.LEFT, padx=10)

        # --- Matplotlibのセットアップ ---
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(self.canvas, self.master, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # --- プロットの初期化 ---
        self.points = self.ax.scatter(self.df['x'], self.df['y'], c='blue', picker=5)
        self.line, = self.ax.plot(self.df['x'], self.df['y'], 'b-', alpha=0.3) # 点の繋がりを可視化
        self.ax.set_title("ウェイポイント (x, y)")
        self.ax.set_xlabel("X座標")
        self.ax.set_ylabel("Y座標")
        self.ax.grid(True)
        self.ax.set_aspect('equal', adjustable='box')

        # --- イベント接続 ---
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('button_press_event', self.on_right_click)

    def on_pick(self, event):
        """点をクリックした際のコールバック。2点を選択して区間を決定する。"""
        if event.artist != self.points:
            return
        
        ind = event.ind[0]
        if ind in self.selection_endpoints: # 同じ点を再度クリックした場合は無視
             return

        self.selection_endpoints.append(ind)
        
        if len(self.selection_endpoints) == 1: # 1点目選択
            self.selected_indices = [ind]
        elif len(self.selection_endpoints) == 2: # 2点目選択
            n = len(self.df)
            p1, p2 = sorted(self.selection_endpoints)
            
            # 経路1（インデックス昇順）
            path1 = list(range(p1, p2 + 1))
            # 経路2（インデックスがループする逆順）
            path2 = list(range(p2, n)) + list(range(0, p1 + 1))
            
            # 短い方の経路を選択
            if len(path1) <= len(path2):
                self.selected_indices = path1
            else:
                self.selected_indices = path2
            
            self.selection_endpoints.clear() # 次の選択のためにリセット

        self.update_plot()

    def on_right_click(self, event):
        """右クリックで選択した点をカーソル位置へ移動させる。"""
        if event.button != 3 or event.inaxes != self.ax or not self.selected_indices:
            return

        target_x, target_y = event.xdata, event.ydata
        move_factor = self.move_factor.get()

        current_coords = self.points.get_offsets()[self.selected_indices]
        new_coords = current_coords + (np.array([target_x, target_y]) - current_coords) * move_factor

        self.df.loc[self.selected_indices, 'x'] = new_coords[:, 0]
        self.df.loc[self.selected_indices, 'y'] = new_coords[:, 1]
        
        self.update_plot()

    def smooth_selection(self):
        """選択された点とその周辺の点に平滑化を適用する。"""
        if not self.selected_indices:
            messagebox.showwarning("警告", "平滑化する点が選択されていません。")
            return

        n = len(self.df)
        padding = self.smooth_padding.get()
        window_half_size = self.smooth_window.get()
        
        # 平滑化を適用するインデックスのセットを作成（選択点＋周辺点）
        indices_to_smooth = set()
        for idx in self.selected_indices:
            for i in range(-padding, padding + 1):
                indices_to_smooth.add((idx + i) % n)

        # 元のデータコピーに対して計算を行う
        original_df = self.df.copy()
        
        for i in sorted(list(indices_to_smooth)):
            # 移動平均を計算するためのインデックス窓（円環を考慮）
            window_indices = [(i + j) % n for j in range(-window_half_size, window_half_size + 1)]
            
            # 新しい座標を計算
            new_x = original_df.loc[window_indices, 'x'].mean()
            new_y = original_df.loc[window_indices, 'y'].mean()
            
            # DataFrameを更新
            self.df.loc[i, 'x'] = new_x
            self.df.loc[i, 'y'] = new_y

        self.update_plot()
        messagebox.showinfo("成功", f"{len(indices_to_smooth)}個の点を平滑化しました。")

    def update_plot(self):
        """プロットの表示を現在のデータに合わせて更新する。"""
        # 点の色を更新
        colors = np.array(['blue'] * len(self.df))
        if self.selected_indices:
            colors[self.selected_indices] = 'red'
        self.points.set_color(colors)
        
        # 点の座標を更新
        all_offsets = self.df[['x', 'y']].to_numpy()
        self.points.set_offsets(all_offsets)
        
        # 線の座標を更新
        line_x = np.append(self.df['x'], self.df['x'][0]) # ループを閉じる
        line_y = np.append(self.df['y'], self.df['y'][0])
        self.line.set_data(line_x, line_y)

        self.canvas.draw_idle()

    def save_csv(self):
        """変更をCSVファイルに保存する。"""
        try:
            self.df.to_csv(self.filepath, index=False)
            messagebox.showinfo("成功", f"ウェイポイントを正常に保存しました:\n{self.filepath}")
        except Exception as e:
            messagebox.showerror("エラー", f"ファイルの保存に失敗しました: {e}")

def main():
    # ここでCSVファイルのパスを直接指定します
    # スクリプトと同じディレクトリに 'waypoints.csv' があることを想定
    csv_path = "output.csv"

    root = tk.Tk()
    app = WaypointEditor(root, csv_path)
    root.mainloop()

if __name__ == "__main__":
    main()