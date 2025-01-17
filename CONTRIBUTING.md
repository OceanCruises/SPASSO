# Contributing to SPASSO

:+1::tada: Thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines (not rule) for contributing to the developpement of SPASSO
Feel free to open a new issue when encountering a bug and propose changes to this 
document in a pull request.

#### Table Of Contents

[Code of Conduct](#code-of-conduct)

[What should I know before I get started?](#what-should-i-know-before-i-get-started)

[How Can I Contribute?](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Your First Code Contribution](#your-first-code-contribution)

## Code of Conduct

SPASSO is a free software designed to support scientist in planning their oceanographic
cruise. This project benefited from the inputs of many contributors over the years 
and the authors wish to maintain this vision. Everyone is free to participate to uphold this code
as long as they adopt respectful behaviors. We as contributors and maintainers pledge 
to make participation in our project and our community a harassment-free experience for everyone,
 regardless of age, body size, disability, ethnicity, gender identity and expression,
  level of experience, nationality, personal appearance, race, religion, or sexual identity 
  and orientation.

## What should I know before I get started?

The SPASSO software have been developped and made freely available by researchers 
from their own field experiments. It is thus not designed to fulfill every oceanographer needs.
Please recall that the mainteners are researchers whose job is not exclusively dedicated to this
software developpement. In this context, we as mainteners are greatly thankful for anyone wishing
to contribute to the project by reporting bugs or suggesting enhancements.

## How Can I Contribute?

### Reporting Bugs

Following these guidelines helps maintainers and the community understand your report :pencil:, reproduce the behavior :computer: :computer:, and find related reports :mag_right:.

Before creating bug reports, please check in [Issues](https://github.com/OceanCruises/SPASSO/issues) as you might find out that you don't need to create one. When you are creating a bug report, please [include as many details as possible](#how-do-i-submit-a-good-bug-report).

> **Note:** If you find a **Closed** issue that seems like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one.

#### How Do I Submit A Bug Report?

When you identify a bug, create an issue and provide as much information as you can on the issue and include additional details to help maintainers reproduce the problem:

* **Use a clear and descriptive title** for the issue to identify the problem.
* **Describe the exact steps which reproduce the problem** in as many details as possible. For example, provide your config_ini file.
* **Describe the behavior you observed after following the steps**. Include screenshots with the errors occuring and/or image outputs.

Provide more context by answering these questions:

* **Did the problem start happening recently** or was this always a problem? This configuration used to work but not anymore ?
* **Which version of Python are you using?** 
* **What's the name and version of the OS you're using**?
...

### Suggesting Enhancements

Everyone is welcomed to suggest including completely new features and minor improvements to existing functionality.

Please include as many detail as possible to describe the modifications/new option. 

Enhancement suggestions have to be published in Issues and must include at least the following information:

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Explain why this enhancement would be useful** to most SPASSO users and/or specific oceanographers community.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**. Include copy/pasteable snippets which you use in those examples
* **Specify the name and version of the OS and Python you're using.**

### Your First Code Contribution
The SPASSO code can be developped locally. Modifications (if relevant for a wider community) can be included in SPASSO. Developpers need to create a Pull Request describing in detail the modifications included in their own version. Please create a config Example (similar to WMedSeaExample) so that maintainers can test before accepting or declining the pull request.
The reviewer(s) may ask you to complete additional design work, tests, or other changes before your pull request can be ultimately accepted.
This process will help maintaining SPASSO's quality, ability to run properly on different machine and fix problems that are important to users.


#### Issues and Pull Request Labels
In the following tables, please find a list of labels that contributors might use to help tracking and managing issues and pull requests.

| Label name | Description |
| --- | --- |
| `enhancement` | Feature requests. |
| `bug` | Confirmed bugs or reports that are very likely to be bugs. |
| `question` | Questions more than bug reports or feature requests (e.g. how do I do X). |
| `feedback` | General feedback more than bug reports or feature requests. |
| `help-wanted` | Users would appreciate help from the community in resolving these issues. |
| `option-idea` | Feature request which might be good candidates for new option to be added to SPASSO. |

| Topic category |  Description |
| --- | --- |
| `download` | Related to data download. |
| `documentation` | Related to any type of documentation. |
| `lagrangian` | Related to Lagrangian diagnostic computation. |
| `eulerian` | Related to Eulerian diagnostic computation. |
| `trajectories` | Related to numerical particle trajectory computation. |
| `plot` | Related to data plotting. |