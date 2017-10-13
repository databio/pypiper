# Design decision notes

## Terms
- **Stage** or **phase**: an arbitrarily defined logical processing *step* or 
*unit* of operation(s) within a pipeline (e.g., read trimming or peak calling
- **Checkpoint**: closely tied to the notion of a *stage* or *phase*, a 
checkpoint represents a point in a pipeline that the author has deemed 
as sufficiently significant to warrant designation.


## Classes

### `Pipeline`

Since a pipeline author detemines how to compose logical units, processing 
steps, or phases to define the notion of a pipeline, this class is 
inherently abstract. We prefer to be able to impose and enforce the requirement 
for stage definitions up front. This precludes the definition or creation of 
a `Pipeline` without stages as we declare `stages` as an `abc.abstractproperty` 
in the definition of `Pipeline`. This also permits us to validate the 
stage definitions up front, at time of pipeline creation rather than waiting 
until invocation of something like `run`. A further benefit of this design is 
the ability to store the parsed, validated form of the stage definitions 
obtained during instance construction. This eliminates a potential need 
to pass the stage definitions among methods for which they're needed, thereby 
simplifying our function signatures.

Not only do we want to provide a simple framework in 
which processing stage/phases may be enumerated and defined in sequence, 
but we also want to facilitate allow the possibility of non-sequential stages 
to be defined by the pipeline author. In the context of say, testing multiple 
alternative ways to do the same conceptual task (e.g., read trimming or peak 
calling) within the same pipeline, in early pipeline development, it's 
particularly likely that the desire to define unordered stages may arise.

Additionally, it would be nice to support varying degrees of expressive power 
and simplicity. To some extent, this is likely to present a trade-off, with 
greater expressive power coming at the expense of implementation simplicity 
for a developer who wishes to implement/extend `Pipeline`. Possibilities 
for some of the "levels" of simplicity and power include but are not limited to:

### `Stage`


## Checkpointing complexity

### Direct pipeline file writes
- In the most basic case, the pipeline may d
