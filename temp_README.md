# [{{package}}]({{url}})

A {{genice}} plugin to show 3D structure models in Jupyter.

version {{version}}

## Requirements

{% for i in requires %}
* {{i}}
{%- endfor %}

## Installation from PyPI

```shell
% pip install git+{{url}}
```

## Manual Installation

### System-wide installation

```shell
% make install
```

### Private installation

Copy the files in {{base}}/formats/ into your local formats/ folder.

## Usage

Open jupyter.ipynb in the Jupyter notebook.

{%- filter indent %}
    {{usage}}
{%- endfilter %}



<!-- ## Test in place

```shell
% make test
``` -->
