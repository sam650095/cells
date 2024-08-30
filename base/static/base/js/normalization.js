// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  await normalization();
  await merged();
  document.getElementById("result").classList.remove("hidden");
  // next page btn
  toggleLoading(false, "processbutton");
  document.getElementById("nextbtn").classList.remove("hidden");
}
async function normalization() {
  const csrftoken = getCookie("csrftoken");
  const form = document.getElementById("normalform");
  const formData = new FormData(form);
  const process_result = await fetchAPI("/api/normal", formData, csrftoken);
  console.log(process_result);
  // addata_results div
  let adata_results = document.getElementById("n_adata_results");
  adata_results.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add("space-y-1", "text-gray-500", "list-disc", "list-inside");
  process_result["adata_results"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  adata_results.appendChild(ul);
}
async function merged() {
  const csrftoken = getCookie("csrftoken");
  const process_result = await fetchAPI("/api/merge", 0, csrftoken);
  console.log(process_result);
  // addata_results div
  let adata_results = document.getElementById("m_adata_results");
  adata_results.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add("space-y-1", "text-gray-500", "list-disc", "list-inside");
  process_result["adata_results"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  adata_results.appendChild(ul);
}
