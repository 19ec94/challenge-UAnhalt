# 01_combinatorics: Distinct Colourings of a Cube

This repository provides code to compute the number of **distinct ways to colour the six faces of a cube** using up to `c` colours, with each colour used at most `m` times. Colourings that are the same up to cube rotations are counted only once, using **Burnside’s Lemma** and the **cycle index polynomial**.

## Features

- Handles up to `c` colours, each appearing at most `m` times (with `1 ≤ m ≤ 6`).
- Counts only genuinely distinct colourings under cube rotation (removing duplicates due to symmetrical arrangements).
- Shows all valid distributions of colours and their associated counts.

## Requirements

- Python 3
- sympy library

### Install requirements

If you have a `requirements.txt` file, run:

```bash
pip install -r requirements.txt
```
Or, just install sympy directly:

```bash
pip install sympy
```

## Usage

Run the script with:

```bash
python3 prob1.py
```

You will be prompted for:

- The number of colours (`c`) to use (integer ≥ 1).
- The maximum number of times any single colour may be used on the cube’s faces (`m`, with `1 ≤ m ≤ 6`).

**Example session:**
```bash
Enter the number of colours (at least 1) to paint the cube's facets:
Number of colours (c) = 3

Enter the maximum number of repetitions allowed per colour (1 to 6):
Maximum repetitions per colour (m) = 3

Using the colour combination (c1:3, c2:3, c3:0), there are 2 distinct ways to colour the cube.

Using the colour combination (c1:3, c2:2, c3:1), there are 3 distinct ways to colour the cube.

...

Total number of distinct colourings with c=3 colours and max repetitions per colour m=3 is: 30 
```


## How does this work?

- Computes the **cycle index polynomial** for the rotational symmetries of the cube.
- Uses symbolic algebra (`sympy`) to expand and evaluate this polynomial subject to your colouring constraints.
- Returns for each valid distribution the number of distinct (non-rotationally-equivalent) colourings.

