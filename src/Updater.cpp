#include <IMP/domino3/Updater.h>

IMPDOMINO3_BEGIN_NAMESPACE

Updater::Updater(const NodesTemp &nodes,
                 std::string name):
  base::Object(name), nodes_(nodes) {
}

void Updater::update(unsigned int iterations) {
  for (unsigned int i = 0; i < iterations; ++i) {
    do_update();
  }
}

IMPDOMINO3_END_NAMESPACE
