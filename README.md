# unfolder: a software to unfold a mesh

## Objective
You goal is to provide an implementation of a couple of unfolding heuristics from 'Schlickenrieder, Wolfram. "Nets of polyhedra." Master's Thesis, Technische Universität Berlin (1997).'
* You can download his thesis from [this link](https://scholar.google.com/scholar?q=Nets%20of%20polyhedra)

## Compile and run
* To compile and run the provided code, please consult the details in the [wiki page](https://github.com/jmlien/unfolder/wiki)

## What to do?
* Understand the code. In particular the files: "Splitter.h/cpp". You will need to provide two functions in these files. Search for "TODO" near the bottom of the files. There are the only files that you have to modify for this assignment.
* Understand "Steepest edge cut tree" and "Flat edge unfolding" in Schlickenrieder's work. Find their corresponding implementation in Splitter.h/cpp and understand the implementation. 
* Implement two heuristics in the classes called MySplitter01 and MySplitter02 in Splitter.h/cpp
  * these two heuristics cannot be one of those that are implemented in Splitter.h/cpp
  * the code is setup in the way that it would be easiest for you to implement the heuristics in "4.7 Direction unfolding" in Schlickenrieder's thesis
* Test your implementation on **5 convex** models and **5 non-convex** models. Do the same using  "Steepest edge cut tree" and "Flat edge unfolding". Report your foundings.

## How to submit?
* Send me an url pointing to your github repository that contains your implementation and the report.
  * Put your report in "unfolder/your-report" directory.
  * You can use your favorite text editor to make the report (MS word, latex, notepad, OSX page, google doc, etc)
* If you are not familiar with github, please consult this [tutorial](https://guides.github.com/activities/hello-world/)
