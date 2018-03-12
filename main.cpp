#include "map_mainwindow.h"
#include <QApplication>

#include "C3.h"
#include "utils.h"
/*! \namespace std
*
* \brief Espace de nommage standard
*
 * Utilise le namespace de la bibliothèque standard
 */
using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MAP_MainWindow w;
    w.show();

   return a.exec();
}
