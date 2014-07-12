This module is for experimenting with a more traditional inferential approach. The computational model is:
- a graph with a node for each piece of data, sharing edges when they share variables
- marginals for each variable, with a copy at each node that uses it
- message passing between nodes when they update the marginals

Note that due to the heavy use of SSE intrinsics, this module currently only
works with gcc.

# Info

_Author(s)_: Daniel Russel, Martin Steinegger

_Maintainer_: `benmwebb`

_License_: None

_Publications_:
- None
