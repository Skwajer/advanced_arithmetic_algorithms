#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include <QVector>

class ExponentPair {
public:
    int alpha, beta;
};

class MonomialIdeal {
public:
    void setGenerators(const QVector<ExponentPair>& gens);
    bool contains(int m, int n) const;
    bool canBeInRemainder(int m, int n) const;
private:
    QVector<ExponentPair> generators;
};

QT_BEGIN_NAMESPACE
class QTableWidget;
class QLineEdit;
class QPushButton;
class QRadioButton;
class QGraphicsView;
class QGraphicsScene;
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void addGenerator();
    void removeSelected();
    void clearAll();
    void updateVisualization();

private:
    void setupUI();

    QTableWidget *table;
    QLineEdit *editAlpha;
    QLineEdit *editBeta;
    QPushButton *btnAdd;
    QPushButton *btnRemove;
    QPushButton *btnClear;
    QPushButton *btnUpdate;
    QRadioButton *radioIdeal;
    QRadioButton *radioRemainder;
    QGraphicsView *view;
    QGraphicsScene *scene;

    MonomialIdeal ideal;
};

#endif
