---
layout: default
title: software
---

### trajectory alignment software

[Download](https://github.com/apicco/trajectory_alignment/archive/master.zip)

[github repository](https://github.com/apicco/trajectory_alignment/)

{% for p in site.pages %}
	{% if p.menu == "wiki" %}
		{% if p.wikiPageName == "Home" %}
        <li><a class="post-link" href="{{ p.url | prepend: site.baseurl }}"> Documentation </a></li>
        {% endif %}
    {% endif %}
{% endfor %}

References:

* [Picco, _et al._, 2015](http://dx.doi.org/10.7554/eLife.04535)

* [Kukulski, Picco, _et al._](http://dx.doi.org/10.7554/eLife.16036)
