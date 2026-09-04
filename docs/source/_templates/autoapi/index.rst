.. This file is the AutoAPI landing page template. It overrides the one shipped
   with sphinx-autoapi (see ``autoapi_template_dir`` in conf.py). Edit the title
   and the prose below freely. Only the Jinja loops are machinery, everything
   else is plain reStructuredText. Note that Jinja delimiters are parsed even
   inside reStructuredText comments such as this one, so do not write any here.

Feature measurements
====================

Each module of the ``features`` folder defines a single class, dedicated to the measurement of one specific feature.
All these classes are called by the :class:`Measure <measure.Measure>` class, which orchestrates the featurisation of a collection of protein structures.
The table below lists the feature classes currently available.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Class
     - Description
{% for page in pages|selectattr("is_top_level_object")|sort(attribute="name") %}
{% for klass in page.classes %}
   * - :py:class:`{{ klass.id }}`
     - {{ klass.docstring|first_sentence or "*(undocumented)*" }}
{% endfor %}
{% if not page.classes %}
   * - :doc:`{{ page.include_path }}`
     - {{ page.docstring|first_sentence or "*(undocumented)*" }}
{% endif %}
{% endfor %}

.. toctree::
   :titlesonly:
   :hidden:

{% for page in pages|selectattr("is_top_level_object")|sort(attribute="name") %}
   {{ page.include_path }}
{% endfor %}
