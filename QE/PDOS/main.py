import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QPushButton,
    QFileDialog, QLabel, QLineEdit, QVBoxLayout, QHBoxLayout,
    QListWidget, QMessageBox, QInputDialog
)
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


def read_scf_fermi(path):
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if "fermi" in line.lower():
                tokens = line.replace("=", " ").split()
                for t in tokens:
                    try:
                        return float(t)
                    except:
                        pass
    raise ValueError("Fermi energy not found in .out file.")


class PDOSApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PDOS Plotter (Editable Legend)")
        self.resize(600, 700)

        self.pdos_files = []
        self.labels = []  
        self.Ef = 0.0
        self.figure = None
        self.ax = None

        self.init_ui()
        self.init_plot_window()

    def init_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        self.btn_add_file = QPushButton("Add PDOS File(s)")
        self.btn_remove_file = QPushButton("Remove Selected File")
        self.btn_load_fermi = QPushButton("Load .out & Shift Fermi")
        self.btn_plot = QPushButton("Plot / Update")
        self.btn_save = QPushButton("Save Figure")

        self.btn_add_file.clicked.connect(self.add_pdos_files)
        self.btn_remove_file.clicked.connect(self.remove_selected_file)
        self.btn_load_fermi.clicked.connect(self.load_fermi_and_shift)
        self.btn_plot.clicked.connect(self.plot)
        self.btn_save.clicked.connect(self.save_figure)
        self.file_list = QListWidget()
        self.title_edit = QLineEdit("PDOS Plot")
        self.x_label = QLineEdit("Energy (E - Ef) [eV]")
        self.y_label = QLineEdit("PDOS / DOS (a.u.)")
        self.xmin = QLineEdit("-5")
        self.xmax = QLineEdit("5")
        self.ymin = QLineEdit("-30")
        self.ymax = QLineEdit("30")
        self.fig_w = QLineEdit("7")
        self.fig_h = QLineEdit("6")

        # Layout
        control = QVBoxLayout()
        control.addWidget(self.btn_add_file)
        control.addWidget(self.btn_remove_file)
        control.addWidget(self.btn_load_fermi)
        control.addWidget(QLabel("Loaded PDOS Files:"))
        control.addWidget(self.file_list)

        control.addWidget(QLabel("Title"))
        control.addWidget(self.title_edit)
        control.addWidget(QLabel("X Label"))
        control.addWidget(self.x_label)
        control.addWidget(QLabel("Y Label"))
        control.addWidget(self.y_label)

        range_layout = QHBoxLayout()
        for lbl, widget in zip(["Xmin", "Xmax", "Ymin", "Ymax"],
                               [self.xmin, self.xmax, self.ymin, self.ymax]):
            range_layout.addWidget(QLabel(lbl))
            range_layout.addWidget(widget)
        control.addLayout(range_layout)

        size_layout = QHBoxLayout()
        size_layout.addWidget(QLabel("Fig Width"))
        size_layout.addWidget(self.fig_w)
        size_layout.addWidget(QLabel("Fig Height"))
        size_layout.addWidget(self.fig_h)
        control.addLayout(size_layout)

        control.addWidget(self.btn_plot)
        control.addWidget(self.btn_save)

        central.setLayout(control)

    def init_plot_window(self):
        self.figure, self.ax = plt.subplots(
            figsize=(float(self.fig_w.text()), float(self.fig_h.text()))
        )
        self.figure.canvas.manager.set_window_title("PDOS Plot")
        self.figure.canvas.mpl_connect("pick_event", self.on_legend_pick)
        plt.ion()
        self.figure.show()

    def add_pdos_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Select PDOS files", "", "Data Files (*.dat *.txt)"
        )
        for p in paths:
            if p not in self.pdos_files:
                self.pdos_files.append(p)
                self.labels.append(p.split("/")[-1])
                self.file_list.addItem(p.split("/")[-1])

    def remove_selected_file(self):
        row = self.file_list.currentRow()
        if row >= 0:
            self.pdos_files.pop(row)
            self.labels.pop(row)
            self.file_list.takeItem(row)
        else:
            QMessageBox.warning(self, "Warning", "Select a file to remove.")

    def load_fermi_and_shift(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select SCF .out file", "", "Output Files (*.out *.txt)"
        )
        if not path:
            return
        try:
            self.Ef = read_scf_fermi(path)
            QMessageBox.information(
                self, "Fermi Energy Loaded",
                f"Ef = {self.Ef:.4f} eV\n(All curves shifted)"
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def plot(self):
        if not self.pdos_files:
            QMessageBox.warning(self, "Warning", "No PDOS files loaded.")
            return

        self.ax.clear()
        self.figure.set_size_inches(float(self.fig_w.text()),
                                    float(self.fig_h.text()))

        self.lines = []
        for path, label in zip(self.pdos_files, self.labels):
            try:
                data = np.loadtxt(path, skiprows=1)
                energy = data[:, 0] - self.Ef
                up = data[:, 1]
                down = -data[:, 2]

                line_up, = self.ax.plot(energy, up, label=label, picker=True)
                self.ax.plot(energy, down)
                self.ax.fill_between(energy, 0, up, alpha=0.25)
                self.ax.fill_between(energy, 0, down, alpha=0.25)

                self.lines.append(line_up)
            except Exception as e:
                print(f"Failed loading {path}: {e}")

        self.ax.axvline(0, linestyle="--", linewidth=1)
        self.ax.set_title(self.title_edit.text())
        self.ax.set_xlabel(self.x_label.text())
        self.ax.set_ylabel(self.y_label.text())
        self.ax.set_xlim(float(self.xmin.text()), float(self.xmax.text()))
        self.ax.set_ylim(float(self.ymin.text()), float(self.ymax.text()))
        self.leg = self.ax.legend(fontsize=8)
        self.leg.set_draggable(True)
        self.ax.grid(alpha=0.3)
        self.figure.canvas.draw_idle()

    def on_legend_pick(self, event):
        if not isinstance(event.artist, plt.Text):
            return
        text_obj = event.artist
        idx = None
        for i, line in enumerate(self.lines):
            if text_obj.get_text() == line.get_label():
                idx = i
                break
        if idx is None:
            return
        new_label, ok = QInputDialog.getText(self, "Edit Legend Label",
                                             "New label:", text=self.labels[idx])
        if ok and new_label:
            self.labels[idx] = new_label
            self.plot() 

    def save_figure(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", "pdos.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        if path:
            self.figure.savefig(path, dpi=400, bbox_inches="tight")
            QMessageBox.information(self, "Saved", f"Saved to:\n{path}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = PDOSApp()
    win.show()
    sys.exit(app.exec_())
