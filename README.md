This module is for experimenting with a more traditional inferential approach. The computational model is:
- a graph with a node for each piece of data, sharing edges when they share variables
- marginals for each variable, with a copy at each node that uses it
- message passing between nodes when they update the marginals

Note that due to the heavy use of SSE intrinsics, this module currently only
works with gcc.

# Info

_Author(s)_: Daniel Russel, Martin Steinegger

_Maintainer_: `benmwebb`

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
- None
