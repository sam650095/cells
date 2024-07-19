// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  const csrftoken = getCookie("csrftoken");
  const process_result = await fetchAPI("/api/merge", 0, csrftoken);
  console.log(process_result);
  // addata_results div
  let adata_results = document.getElementById("adata_results");
  adata_results.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add("space-y-1", "text-gray-500", "list-disc", "list-inside");
  process_result["adata_results"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  adata_results.appendChild(ul);
  // next page btn
  toggleLoading(false, "processbutton");
  document.getElementById("nextbtn").classList.remove("hidden");
}
