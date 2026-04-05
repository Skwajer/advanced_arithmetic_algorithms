#include "mainwindow.hpp"
#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QLineEdit>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QRadioButton>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsTextItem>
#include <QHeaderView>
#include <QMessageBox>
#include <QFrame>

void MonomialIdeal::setGenerators(const QVector<ExponentPair>& gens) {
    generators = gens;
}

bool MonomialIdeal::contains(int m, int n) const {
    for (const auto& g : generators) {
        if (m >= g.alpha && n >= g.beta) return true;
    }
    return false;
}

bool MonomialIdeal::canBeInRemainder(int m, int n) const {
    for (const auto& g : generators) {
        if (m >= g.alpha && n >= g.beta) return false;
    }
    return true;
}

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    setupUI();
    updateVisualization();
}

MainWindow::~MainWindow() {}

void MainWindow::setupUI() {
    setWindowTitle("Мономиальный идеал");
    resize(900, 700);

    QWidget *central = new QWidget(this);
    setCentralWidget(central);

    QHBoxLayout *mainLayout = new QHBoxLayout(central);

    QWidget *leftPanel = new QWidget();
    QVBoxLayout *leftLayout = new QVBoxLayout(leftPanel);
    leftLayout->setAlignment(Qt::AlignTop);

    leftLayout->addWidget(new QLabel("Образующие (α, β):"));

    table = new QTableWidget(0, 2);
    QStringList headers;
    headers << "α" << "β";
    table->setHorizontalHeaderLabels(headers);
    table->horizontalHeader()->setStretchLastSection(true);
    leftLayout->addWidget(table);

    QHBoxLayout *addRow = new QHBoxLayout();
    editAlpha = new QLineEdit();
    editAlpha->setPlaceholderText("α");
    editBeta = new QLineEdit();
    editBeta->setPlaceholderText("β");
    btnAdd = new QPushButton("+");
    addRow->addWidget(editAlpha);
    addRow->addWidget(editBeta);
    addRow->addWidget(btnAdd);
    leftLayout->addLayout(addRow);

    btnRemove = new QPushButton("Удалить выбранную");
    leftLayout->addWidget(btnRemove);

    btnClear = new QPushButton("Очистить всё");
    leftLayout->addWidget(btnClear);

    leftLayout->addWidget(new QFrame());

    leftLayout->addWidget(new QLabel("Режим:"));
    radioIdeal = new QRadioButton("Мономы в I");
    radioRemainder = new QRadioButton("Мономы в остатке");
    radioIdeal->setChecked(true);
    leftLayout->addWidget(radioIdeal);
    leftLayout->addWidget(radioRemainder);

    btnUpdate = new QPushButton("Обновить");
    leftLayout->addWidget(btnUpdate);

    leftLayout->addStretch();

    view = new QGraphicsView();
    scene = new QGraphicsScene(this);
    view->setScene(scene);
    view->setRenderHint(QPainter::Antialiasing);

    mainLayout->addWidget(leftPanel, 1);
    mainLayout->addWidget(view, 2);

    connect(btnAdd, &QPushButton::clicked, this, &MainWindow::addGenerator);
    connect(btnRemove, &QPushButton::clicked, this, &MainWindow::removeSelected);
    connect(btnClear, &QPushButton::clicked, this, &MainWindow::clearAll);
    connect(btnUpdate, &QPushButton::clicked, this, &MainWindow::updateVisualization);
}

void MainWindow::addGenerator() {
    bool ok1, ok2;
    int a = editAlpha->text().toInt(&ok1);
    int b = editBeta->text().toInt(&ok2);

    if (!ok1 || !ok2 || a < 0 || b < 0) {
        QMessageBox::warning(this, "Ошибка", "Введите целые числа ≥ 0");
        return;
    }

    int row = table->rowCount();
    table->insertRow(row);
    table->setItem(row, 0, new QTableWidgetItem(QString::number(a)));
    table->setItem(row, 1, new QTableWidgetItem(QString::number(b)));

    editAlpha->clear();
    editBeta->clear();
    updateVisualization();
}

void MainWindow::removeSelected() {
    int row = table->currentRow();
    if (row >= 0) {
        table->removeRow(row);
        updateVisualization();
    }
}

void MainWindow::clearAll() {
    table->setRowCount(0);
    updateVisualization();
}

void MainWindow::updateVisualization() {
    QVector<ExponentPair> gens;
    for (int r = 0; r < table->rowCount(); ++r) {
        QTableWidgetItem *ia = table->item(r, 0);
        QTableWidgetItem *ib = table->item(r, 1);
        if (ia && ib) {
            ExponentPair p;
            p.alpha = ia->text().toInt();
            p.beta = ib->text().toInt();
            gens.append(p);
        }
    }
    ideal.setGenerators(gens);

    const int N = 20;
    const int cell = 20;
    const int offX = 40;
    const int offY = 40;

    scene->clear();

    QPen pen(Qt::gray, 0.5);

    QBrush brushIdeal(QColor(180, 210, 240));
    QBrush brushRemainder(QColor(200, 230, 190));
    QBrush brushWhite(Qt::white);

    bool modeIdeal = radioIdeal->isChecked();

    for (int m = 0; m <= N; ++m) {
        for (int n = 0; n <= N; ++n) {
            int x = offX + m * cell;
            int y = offY + (N - n) * cell;

            bool inI = ideal.contains(m, n);
            bool inRem = ideal.canBeInRemainder(m, n);

            QBrush brush;
            if (modeIdeal) {
                brush = inI ? brushIdeal : brushWhite;
            } else {
                brush = inRem ? brushRemainder : brushWhite;
            }

            scene->addRect(x, y, cell-1, cell-1, pen, brush);
        }
    }

    for (const auto& g : gens) {
        int m = g.alpha;
        int n = g.beta;
        if (m >= 0 && m <= N && n >= 0 && n <= N) {
            int x = offX + m * cell;
            int y = offY + (N - n) * cell;

            QPen redPen(QColor(200, 50, 50), 2);
            scene->addRect(x, y, cell-1, cell-1, redPen);

            int cx = x + cell/2;
            int cy = y + cell/2;
            scene->addEllipse(cx - 5, cy - 5, 10, 10, QPen(Qt::black), QBrush(QColor(50, 80, 150)));
        }
    }

    QFont font("Arial", 8);
    for (int i = 0; i <= N; ++i) {
        QGraphicsTextItem *tx = scene->addText(QString::number(i), font);
        tx->setPos(offX + i * cell, offY - 15);

        QGraphicsTextItem *ty = scene->addText(QString::number(i), font);
        ty->setPos(offX - 20, offY + (N - i) * cell);
    }

    QGraphicsTextItem *lx = scene->addText("m (степень x)", font);
    lx->setPos(offX + N * cell / 2, offY - 30);

    QGraphicsTextItem *ly = scene->addText("n (степень y)", font);
    ly->setPos(offX - 35, offY + N * cell + 15);
}
