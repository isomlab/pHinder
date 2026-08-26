"""A budget for the jostle loops that move vertices into general position.

Every one of those loops was written as `while not <test>: v.jostle()`, with no
limit. That shape is a silent hang: if the test cannot be satisfied by moving
the vertex the loop moves, the process spins at full CPU forever, writing
nothing and reporting nothing.

That was not hypothetical. gp2D compared only x and y, so two vertices sharing
those two coordinates made it report collinearity for every possible third
vertex; the loop jostled the third vertex, which could not affect the other two,
and ran to millions of iterations while the vertex wandered hundreds of angstroms
from the structure. The predicate is fixed, but the loop shape is worth fixing
too: a bounded loop that names what failed is recoverable, and an unbounded one
is not.

The limit is deliberately far above what any real structure needs -- ordinary
inputs clear these tests in a handful of jostles, or none.
"""

MAX_JOSTLES = 1000


class GeneralPositionError(RuntimeError):
    """Jostling a vertex could not bring a simplex into general position."""


class JostleBudget:
    """Counts the jostles spent by one general-position loop.

    Use it in place of a bare ``v.jostle()`` inside such a loop::

        budget = JostleBudget("v1, v2 and v3 are collinear")
        while not vGp:
            vGp = gp2D(v1, v2, v3)
            if not vGp:
                budget.jostle(v3)

    Keyword arguments are passed through to the vertex's own jostle, so
    ``budget.jostle(v3, parabaloid=True)`` works as before.
    """

    def __init__(self, description, max_jostles=MAX_JOSTLES):
        self.description = description
        self.max_jostles = max_jostles
        self.attempts = 0

    def jostle(self, vertex, **kwargs):
        if self.attempts >= self.max_jostles:
            raise GeneralPositionError(
                "%s: still degenerate after %d jostles, so the calculation was "
                "stopped rather than left to spin. The vertex being moved is at "
                "(%.6f, %.6f, %.6f). This means either that the structure is "
                "genuinely degenerate -- many exactly collinear, coplanar or "
                "coincident atoms, which real coordinates rarely are -- or that "
                "the test does not depend on the vertex this loop moves."
                % (self.description, self.max_jostles,
                   vertex.x, vertex.y, vertex.z)
            )
        self.attempts += 1
        vertex.jostle(**kwargs)
