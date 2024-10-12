from functools import partial
from operator import xor
from pathlib import Path
import re
from types import MappingProxyType as MPt
from typing import (
    Any,
    Callable,
    Collection,
    Mapping,
    Optional,
    Sequence,
    TypeVar,
    Union,
)
import warnings

from cytoolz import keyfilter
from lxml import etree
from more_itertools import chunked, nth

BLANKLINE = re.compile(r"\s*\n")
EXTRA_NEWLINES = re.compile(r"(?<=\n) *\n", re.M)
DELETE_TAG = re.compile(r"\s*{delete:(\w.*?)}")
TEMPLATE_TAG = re.compile(r"{.*?}")
XML_COMMENT = re.compile(r"\s*<!--")
XML_TAG = re.compile(r"<.*?>")
STOP_TAG = re.compile(r"\s*{stop:(\w.*?)}")


def _check_deletion(
    line: str,
    delete_tag: Optional[str],
    targetpat: Optional[re.Pattern],
) -> tuple[bool, Optional[str]]:
    # always delete template comments.
    if XML_COMMENT.match(line):
        return False, delete_tag
    if delete_tag is None:  # we are not in deletion mode.
        # case 1: encountered an unmatched stop command. delete it and move on.
        if STOP_TAG.match(line) is not None:
            return False, None
        # case 2: encountered a deletion start tag. check if we should turn
        # deletion mode on (and always delete the deletion command itself)
        if (dmatch := DELETE_TAG.match(line)) is not None:
            if (tmatch := targetpat.match(dmatch.group(1))) is None:
                return False, None
            return False, tmatch.group()
        # case 3: it's just a regular line. write it.
        return True, delete_tag
    # if we reach this block, we are in deletion mode.
    # case 1: did not encounter stop tag. continue deleting.
    if (stopmatch := STOP_TAG.match(line)) is None:
        return False, delete_tag
    # case 2: encountered stop tag. if it matches opening deletion tag,
    # turn off deletion mode. always delete the stop command itself.
    if delete_tag == stopmatch.group(1):
        return False, None
    # case 3: it's just a regular line; delete it.
    return False, delete_tag


def remove_deletion_tags(template):
    return [
        line for line in template
        if DELETE_TAG.match(line) is None and XML_COMMENT.match(line) is None
    ]


def process_deletions(
    lines: Sequence[str], targets: Collection[str]
) -> list[str]:
    output_lines, delete_tag = [], None
    targetpat = re.compile('|'.join(targets))
    for line in lines:
        write_line, delete_tag = _check_deletion(line, delete_tag, targetpat)
        if write_line:
            output_lines.append(line)
    return output_lines


TemplateParameter = TypeVar(
    "TemplateParameter",
    str,
    tuple[str, Union[str, int]],
    tuple[str, Callable[[Any], str]],
    tuple[str, Callable[[Any], str], Union[str, int]],
)


def prep_tagdict(tagdict: Mapping[str, TemplateParameter]) -> dict[str, Any]:
    """
    format a human-written tag dictionary for label conversion.
    legal signatures for values of dictionary:
        key: str -> metaget(key)
        (key: str, func: Callable) -> func(metaget(key))
        (key: str, subkey: Union[str, int])
            -> metaget(key)][subkey] (for quantities etc.)
        (key: str, func: Callable, subkey: Union[str, int)
            -> func(metaget(key)[subkey])
    """
    prepped = {}
    for tag, parameter in tagdict.items():
        if isinstance(parameter, str):
            prepped[tag] = (parameter, None, None)
        elif len(parameter) == 3:
            prepped[tag] = parameter
        elif isinstance(parameter[1], Callable):
            prepped[tag] = (parameter[0], parameter[1], None)
        else:
            prepped[tag] = (parameter[0], None, parameter[1])
    return prepped


def xmlwrap(
    element_name: str, **element_attributes: str
) -> Callable[[str], str]:
    """
    designed for cases in which the parameters present in a source label are
    widely variable even within a single 'product type' and so variability
    cannot be cleanly handled with deletion tags.
    """
    attribute_string = ""
    if len(element_attributes) > 0:
        attribute_string += " "
        attribute_string += " ".join(
            [f"{attr}={value}" for attr, value in element_attributes.items()]
        )
    opening = f"<{element_name}{attribute_string}>"
    closing = f"</{element_name}>"

    def wrap_in_xml(content):
        return f"{opening}{content}{closing}"

    return wrap_in_xml


def skipna(func: Callable[[str, ...], str]) -> Callable[[str, ...], str]:
    """skip not-found or N/A values"""

    def skipper(content, *args, **kwargs):
        if content in ("N/A", None):
            return ""
        return func(content, *args, **kwargs)

    return skipper


def template_tags(template_fn: Union[str, Path]) -> set[str]:
    """return set of tags from a template file"""
    with open(template_fn) as template:
        tags = re.findall(r"{.*?}", template.read())
    return set(tags)


def inside_the_tag(xml_line: str) -> map:
    """
    Get values from inside an XML tag.

    Caution:
        This is not a 'real' general-purpose XML-parsing function! It is
        intended to be fast and easy, but takes no care to be _comprehensive_.
        You may not wish to use it on XML you did not generate yourself.
    """
    return map(partial(nth, n=1), (chunked(XML_TAG.split(xml_line), 3)))


def manipulate_xml_tag(
    label: str, tag: str, manipulator: Callable[[str, str], str]
) -> str:
    """
    Apply `manipulator()` to all values in `label` that appear to be wrapped in
    `tag`.

    Caution:
        This is not a 'real' general-purpose XML-parsing function! It is
        intended to be fast and easy, but takes no care to be _comprehensive_.
        You may not wish to use it on XML you did not generate yourself.
    """
    lines = label.splitlines()
    output_lines = []
    for line in lines:
        if tag in line:
            value = next(inside_the_tag(line))
            output_lines.append(manipulator(line, value))
        else:
            output_lines.append(line)
    return "\n".join(output_lines)


def round_xml_tag(label: str, tag: str, precision: int) -> str:
    def rounder(line: str, value: str) -> str:
        rounded = round(float(value), precision)
        if precision == 0:
            rounded = int(rounded)
        return line.replace(value, str(rounded))

    return manipulate_xml_tag(label, tag, rounder)


from textwrap import wrap
from collections import defaultdict
from typing import Hashable, MutableMapping, MutableSequence


def seqmerge(
    lmap: MutableMapping[Hashable, MutableSequence],
    rmap: MutableMapping[Hashable, MutableSequence]
):
    for k, v in rmap.items():
        lmap[k] = lmap.get(k, []) + v
    return lmap


def sort_by_level(node, level=1):
    nodes = defaultdict(list, {level - 1: [node]})
    for child in node:
        seqmerge(nodes, sort_by_level(child, level + 1))
    return dict(nodes)


def format_processing_instructions(rootnode):
    """etree just does not handle them well."""
    previous = [rootnode.getprevious()]
    while previous[-1] is not None:
        previous.append(previous[-1].getprevious())
    return previous[:-1]


def indent_xml(text, indentation=1, maxwidth=79):
    root = etree.fromstring(text)
    levels = sort_by_level(root)
    istring = ' ' * indentation
    for level, elements in levels.items():
        lstring = istring * level
        for element in elements:
            if (previous := element.getprevious()) is not None:
                previous.tail = f"\n{lstring}"
            elif (parent := element.getparent()) is not None:
                parent.text = parent.text + f"\n{lstring}"
            if element.getnext() is None:
                element.tail = f"\n{istring * (level - 1)}"
            if element.text is None:
                # this condition indicates a nil value
                continue
            text = re.sub(" ?\n( +)?| +", " ", element.text.strip())
            # NOTE: can't wrap LIDS / LIDVIDs
            if len(text) > maxwidth and not text.startswith("urn:"):
                lines = wrap(
                    text,
                    width=maxwidth,
                    initial_indent=f"{lstring}{istring}",
                    subsequent_indent=f"{lstring}{istring}"
                )
                element.text = f"\n{'\n'.join(lines)}\n{lstring}"
            else:
                element.text = text
    instructions = [
        etree.tostring(inst).decode('utf-8').strip()
        for inst in format_processing_instructions(root)
    ]
    body = etree.tostring(root).decode('utf-8')
    return f"{'\n'.join(instructions)}\n{body}"


def replace_nicely_formatted_schema(
    xmlstr: str, output_lines: list[str]
) -> str:
    start, stop = None, None
    for i, l in enumerate(output_lines):
        if start is None and '<Product_' in l:
            start = i
            continue
        if l.strip().startswith('>'):
            stop = i + 1
            break
    if start is None:
        return xmlstr
    if stop is None:
        warnings.warn("Couldn't find nice schema.")
        return xmlstr
    target = output_lines[start].strip()
    nice = re.sub(rf"{target}.*\n", "".join(output_lines[start:stop]), xmlstr)
    if nice == xmlstr:
        warnings.warn("Couldn't find nice schema.")
    return nice


class TagIsNoneError(TypeError):
    pass


def fill_template(
    template_path: Optional[Union[str, Path]] = None,
    template_text: Optional[str] = None,
    deletion_targets: Optional[Sequence[str]] = None,
    associations: Mapping = MPt({}),
    pprint_xml: bool = True,
) -> str:
    if not xor(template_path is None, template_text is None):
        raise TypeError(
            "Exactly one of template_path or template_text must be defined."
        )
    elif template_text is not None:
        template = template_text.splitlines()
    else:
        template = Path(template_path).open().readlines()
    if deletion_targets is not None:
        template = process_deletions(template, deletion_targets)
    else:
        template = remove_deletion_tags(template)
    association_tags = {"{" + k + "}": str(v) for k, v in associations.items()}
    output_lines = []
    for line in template:
        for k, v in keyfilter(lambda k: k in line, association_tags).items():
            if v == "None":
                raise TagIsNoneError(
                    f"refusing to fill tag {k} with None; line is {line}"
                )
            line = line.replace(k, v)
        if BLANKLINE.match(line) is not True:
            output_lines.append(line)
    if pprint_xml is False:
        return EXTRA_NEWLINES.sub("", "\n".join(output_lines))
    try:
        xmlstr = indent_xml("\n".join(output_lines).encode('utf-8'))
        return replace_nicely_formatted_schema(xmlstr, output_lines)
    except (AttributeError, ValueError, KeyError) as ex:
        raise ValueError(
            f"Could not parse XML for pretty-printing\n{type(ex)}: {ex}"
        )